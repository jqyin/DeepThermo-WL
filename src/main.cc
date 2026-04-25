#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>

#include "mpi.h"

#include "alloy.hpp"
#include "allreduce.h"
#include "backend/inference_backend.hpp"
#include "parameter.hpp"
#include "pt.hpp"
#include "rand.hpp"
#include "wanglandau.hpp"

namespace {
int default_gpus_per_node() {
#if defined(DEEPTHERMO_PLATFORM_FRONTIER)
    return 8;
#elif defined(DEEPTHERMO_PLATFORM_PERLMUTTER)
    return 4;
#else
    return 6;
#endif
}
}  // namespace

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiState.nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiState.myrank);
    mpiState.gpus_per_node = default_gpus_per_node();

    if (argc < 3) {
        if (mpiState.myrank == 0) {
            std::fprintf(stderr,
                         "Usage: %s <input-file> <seed>\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int user_seed = std::atoi(argv[2]);
    std::srand(user_seed * (mpiState.myrank + 1));
    shelltimeseed(std::rand() + 19 * mpiState.myrank + 19273);

    try {
        ReadInput(argv[1]);
    } catch (const std::exception& e) {
        if (mpiState.myrank == 0) std::fprintf(stderr, "%s\n", e.what());
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    deepthermo_sim.backend = make_backend();
    try {
        deepthermo_sim.backend->load(alloyState.model_dir,
                                     mpiState.myrank % mpiState.gpus_per_node,
                                     mpiState.myrank, alloyState);
    } catch (const std::exception& e) {
        if (mpiState.myrank == 0) {
            std::fprintf(stderr, "backend (%s) load failed: %s\n",
                         deepthermo_sim.backend->name(), e.what());
        }
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    ini_T(ptState.MTi, ptState.MTf, mpiState.nprocs);
    ini_sys();

    if (mpiState.myrank == 0) ini_alloy(0);
    MPI_Bcast(alloyState.Atom.data(), alloyState.N_3, MPI_SHORT, 0, MPI_COMM_WORLD);
    MPI_Bcast(alloyState.NT.data(), alloyState.NE, MPI_INT, 0, MPI_COMM_WORLD);
    ini_apos();
    alloyState.currEtot = Etot();
    std::printf("Etot = %g\n", alloyState.currEtot);

    // Parallel-tempering warm-up, paper §III-C stage 1.
    parallel_tempering(mpiState.nprocs, ptState.MDROP, ptState.MSAMPS, ptState.MSEP, 0, metropolis);

    double ptEmin = 0.0, ptEmax = 0.0;
    MPI_Allreduce(&alloyState.currEtot, &ptEmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&alloyState.currEtot, &ptEmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (ptEmin / alloyState.N_3 < wlState.WLD1min) wlState.WLD1min = ptEmin / alloyState.N_3;
    if (ptEmax / alloyState.N_3 > wlState.WLD1max) wlState.WLD1max = ptEmax / alloyState.N_3;
    std::printf("ptEmin = %g  ptEmax = %g\n", ptEmin, ptEmax);

    FILE* ofp_run = nullptr;
    FILE* ofp_comm = nullptr;
    std::time_t start = 0;
    if (mpiState.myrank == 0) {
        ofp_run = std::fopen("run.dat", "a");
        ofp_comm = std::fopen("comm.dat", "w");
        start = std::time(nullptr);
        std::fprintf(ofp_run, "#seeds: %d,%d,%d\n", 314159265, 362436069, user_seed);
    }

    wlState.numf = 1;
    const std::size_t msgSize = wlState.D1BINS * sizeof(double);
    const int nsweeps_per_iter = alloyState.N_3 / mpiState.nprocs > 0
                                     ? alloyState.N_3 / mpiState.nprocs
                                     : 1;
    (void)nsweeps_per_iter;  // retained for future tuning

    initWL();
    if (mpiState.myrank == 0) {
        std::printf("Ecut: %f\n",
                    (alloyState.E_0 / wlState.invdWLD1 + wlState.WLD1min) * alloyState.N_3);
    }

    std::chrono::high_resolution_clock::time_point tik, tok;
    double tmp_flat = 0.0;
    for (wlState.lnwlf = wlState.ModFactorInit;
         wlState.lnwlf > wlState.ModFactorFinal;
         wlState.lnwlf = wlState.lnwlf / wlState.IterationFactor) {

        wlState.IterSweeps = 0;
        resetWL();
        tmp_flat = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);

#ifdef GLOBAL_UPDATE
        global_update(alloyState.N_3, wlState.Flatness * wlState.lnwlf);
#endif

        while (tmp_flat <= wlState.Flatness) {
            // θ gate: 1 VAE global update every 2^numf sweeps, once the
            // current energy drops below the shrinking E_0 cutoff.
            const bool theta =
                wlState.TotalSweeps % static_cast<int>(std::pow(2, wlState.numf)) == 0 &&
                alloyState.currEtot <
                    (alloyState.E_0 / wlState.invdWLD1 + wlState.WLD1min) * alloyState.N_3;
            if (theta) vae_update(WLdos);

            sweepWL(alloyState.N_3, WLdos);
            wlState.IterSweeps += 1;
            wlState.TotalSweeps += 1;

            MPI_Barrier(MPI_COMM_WORLD);
            if (wlState.IterSweeps == 10 && mpiState.myrank == 0) {
                tik = std::chrono::high_resolution_clock::now();
            }

#ifdef COMM_RING
            RingAllreduce(wlState.wllngi.data(), msgSize, &wlState.wllngd,
                          mpiState.myrank, mpiState.nprocs);
            RingAllreduce(wlState.wlHi.data(), msgSize, &wlState.wlHd,
                          mpiState.myrank, mpiState.nprocs);
            MPI_Barrier(MPI_COMM_WORLD);
#else
            MPI_Allreduce(wlState.wlHi.data(), wlState.wlHd.data(),
                          wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (wlState.IterSweeps == 10 && mpiState.myrank == 0) {
                tok = std::chrono::high_resolution_clock::now();
            }

            for (int i = 0; i < wlState.D1BINS; ++i) {
                wlState.wlH[i] += wlState.wlHd[i];
                wlState.wllng[i] += wlState.wlHd[i] * wlState.lnwlf;
            }
            std::fill(wlState.wllngi.begin(), wlState.wllngi.end(), 0.0);
            std::fill(wlState.wlHi.begin(), wlState.wlHi.end(), 0);
            tmp_flat = flatWL(WLdos);
            if ((wlState.IterSweeps % 100) == 0 && mpiState.myrank == 0) write_DOS_H();
        }

        if (mpiState.myrank == 0) write_DOS_H();
        if (mpiState.myrank == 0) {
            std::fprintf(ofp_run, "%g\t%d\t%d\t%g\t%g\n",
                         wlState.lnwlf,
                         wlState.TotalSweeps * wlState.D1BINS,
                         wlState.IterSweeps * wlState.D1BINS,
                         tmp_flat,
                         wlState.numbelow_flat);
            std::chrono::duration<double, std::milli> ms_double = tok - tik;
            std::fprintf(ofp_comm,
                         "comm time: %f (ms)\nbandwidth: %f (GB/s)\n",
                         ms_double.count(),
                         2.0 * msgSize / 1e6 / (ms_double.count() > 0 ? ms_double.count() : 1.0));
            std::fflush(ofp_run);
            std::fflush(ofp_comm);
        }
        wlState.numf += 1;
        alloyState.E_0 -= static_cast<int>(0.01 * wlState.D1BINS);
    }

    if (mpiState.myrank == 0) {
        const std::time_t end = std::time(nullptr);
        std::fprintf(ofp_run, "wl time: %ld (s) \n", static_cast<long>(end - start));
        std::fclose(ofp_run);
        std::fclose(ofp_comm);
    }

    // Production run.
    resetWL();
    int count = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    while (count < wlState.ProductionBinSamps) {
        sweepWL(alloyState.N_3, WLproduction);
        MPI_Allreduce(wlState.wlHi.data(), wlState.wlHd.data(),
                      wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < wlState.D1BINS; ++i) {
            wlState.wlH[i] += wlState.wlHd[i];
            wlState.wllng[i] += wlState.wlHd[i] * wlState.lnwlf;
        }
        std::fill(wlState.wlHi.begin(), wlState.wlHi.end(), 0);
        count++;
        wlState.IterSweeps += 1;
        wlState.TotalSweeps += 1;
    }

    std::fill(wlState.wllngi.begin(), wlState.wllngi.end(), 0.0);
    std::fill(wlState.wllngd.begin(), wlState.wllngd.end(), 0.0);
    MPI_Reduce(alloyState.op.data(), wlState.wllngi.data(), wlState.D1BINS,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(alloyState.op2.data(), wlState.wllngd.data(), wlState.D1BINS,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpiState.myrank == 0) {
        write_DOS_H();
        thermoqs();
    }

    MPI_Finalize();
    return 0;
}
