#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "mpi.h"

#include "alloy.hpp"
#include "pt.hpp"
#include "rand.hpp"
#include "wanglandau.hpp"

double flatWL(SamplingMode /*mode*/) {
    int k = 0;
    double min = 1.0e300;
    double avg = 0.0;
    wlState.numbelow_flat = 0.0;

    for (int i = 0; i < wlState.D1BINS; ++i) {
        const double Hi = wlState.wlH[i];
        if (wlState.mask[i] == 1) {
            if (Hi < min) min = Hi;
            if (Hi > 0.0) {
                avg += Hi;
                k += 1;
            }
        }
    }
    avg /= (k > 0 ? 1.0 * k : 1.0);
    if (k == 0) avg = 1.0;

    wlState.numbelow_flat = 0.0;
    return min / avg;
}

void resetWL() {
    std::fill(wlState.wlH.begin(), wlState.wlH.end(), 0.0);
}

void readmask() {
    FILE* ifp = std::fopen("mask.dat", "r");
    if (ifp == nullptr) return;
    for (int k = 0; k < wlState.D1BINS; ++k) {
        int i = 0, tmp = 0;
        std::fscanf(ifp, "%d\t%d\n", &i, &tmp);
        wlState.mask[i] = tmp;
    }
    std::fclose(ifp);
}

void writemask() {
    FILE* ofp = std::fopen("mask.dat", "w");
    for (int i = 0; i < wlState.D1BINS; ++i) {
        std::fprintf(ofp, "%g\t%d\n", i / wlState.invdWLD1 + wlState.WLD1min, wlState.mask[i]);
    }
    std::fclose(ofp);
}

void readg() {
    FILE* ifp = std::fopen("g.dat", "r");
    if (ifp == nullptr) return;

    double tmp = 0.0;
    std::fscanf(ifp, "#%d\t%lg\t%d\t%d\t%lg\n",
                &wlState.numf, &wlState.ModFactorInit,
                &wlState.TotalSweeps, &wlState.IterSweeps, &tmp);
    std::printf("restart from: numf %d  lnwlf %g \n", wlState.numf, wlState.ModFactorInit);
    for (int i = 0; i < wlState.D1BINS; ++i) {
        std::fscanf(ifp, "%lg\t%lg\t%lg\t%lg \n", &tmp, &wlState.wllng[i], &tmp, &tmp);
    }
    std::fclose(ifp);
}

void write_DOS_H() {
    char s1[512];
    std::snprintf(s1, sizeof(s1), "DOS_H_iter%03d.dat", wlState.numf);
    FILE* ofp = std::fopen(s1, "w");

    double maxg = -1.0e300;
    double ming = 1.0e300;
    for (int i = 0; i < wlState.D1BINS; ++i) {
        if ((wlState.wllng[i] > maxg) && (wlState.mask[i] == 1)) maxg = wlState.wllng[i];
        if ((wlState.wllng[i] < ming) && (wlState.mask[i] == 1)) ming = wlState.wllng[i];
    }

    std::fprintf(ofp, "#%d\t%g\t%d\t%d\t%g\n",
                 wlState.numf, wlState.lnwlf,
                 wlState.TotalSweeps, wlState.IterSweeps,
                 flatWL(WLdos));
    for (int i = 0; i < wlState.D1BINS; ++i) {
        std::fprintf(ofp, "%.8f\t%18.10e\t%18.10e\t%g\n",
                     i / wlState.invdWLD1 + wlState.WLD1min,
                     wlState.wllng[i] - maxg,
                     1.0 * wlState.wlH[i],
                     1.0 * wlState.acceptrot[i] / (wlState.attemptrot[i] > 0 ? wlState.attemptrot[i] : 1));
    }
    std::fflush(ofp);
    std::fclose(ofp);
}

void initWL() {
    wlState.TotalSweeps = 1;
    wlState.lnwlf = wlState.ModFactorInit;

    const int bins = static_cast<int>(
        (wlState.WLD1max - wlState.WLD1min) / (wlState.dWLD1 * alloyState.invN));
    wlState.invdWLD1 = 1.0 / (wlState.dWLD1 * alloyState.invN);
    wlState.WLD1max = wlState.WLD1min + bins / wlState.invdWLD1;

    wlState.resize(bins);

    // Order-parameter accumulators live in AlloyState but are bin-indexed,
    // so size them here where the bin count is finalised.
    alloyState.op.assign(bins, 0.0);
    alloyState.op2.assign(bins, 0.0);

    // Relax into the WL energy window with Metropolis sweeps.
    ptState.pT = 100;
    while (alloyState.currEtot * alloyState.invN > wlState.WLD1max) {
        mchybrid(metropolis);
        std::fprintf(stderr, "rank%d Relaxing:  %g of %g \n",
                     mpiState.myrank, alloyState.currEtot * alloyState.invN, wlState.WLD1max);
    }
    if (alloyState.currEtot * alloyState.invN < wlState.WLD1min) {
        std::fprintf(stderr, "currEtot/N: %g  \n", alloyState.currEtot * alloyState.invN);
    }
    wlState.LOWESTE = 0;

    alloyState.E_0 = static_cast<int>(0.2 * wlState.D1BINS);
    readg();
    alloyState.E_0 -= static_cast<int>((wlState.numf - 1) * 0.01 * wlState.D1BINS);

    MPI_Barrier(MPI_COMM_WORLD);
}

void sweepWL(int nsweeps, SamplingMode mode) {
    for (int i = 0; i < nsweeps; ++i) BondSwap(mode);
}

void report_vae(int sweeps, double lnwlf) {
    // Collective: every rank must call this. Counts are cumulative, so
    // differencing successive rows gives the rate within an interval.
    int vloc[2] = {alloyState.attv, alloyState.accv};
    int vsum[2] = {0, 0};
    MPI_Reduce(vloc, vsum, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpiState.myrank != 0) return;

    FILE* fp = std::fopen("vae.dat", "a");
    if (fp == nullptr) return;
    std::fprintf(fp, "%d\t%g\t%d\t%d\t%g\n", sweeps, lnwlf, vsum[0], vsum[1],
                 vsum[0] > 0 ? 1.0 * vsum[1] / vsum[0] : 0.0);
    std::fclose(fp);
}

void global_update(int nsweeps, double Flat) {
    double duration = 0.0;
    double tmp_flat = 0.0;
    int cnt = 0;
    FILE* fp_time = nullptr;
    std::chrono::high_resolution_clock::time_point start, end;

    if (wlState.numf == 2 && mpiState.myrank == 0) {
        fp_time = std::fopen("infer.dat", "w");
    }

    while (tmp_flat <= Flat) {
        sweepWL(nsweeps, WLprior);
        if (wlState.numf == 2) start = std::chrono::high_resolution_clock::now();
        vae_update(WLprior);
        if (wlState.numf == 2) {
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> ms_double = end - start;
            duration += ms_double.count();
            cnt++;
        }

        wlState.IterSweeps += 1;
        wlState.TotalSweeps += 1;
        if ((wlState.IterSweeps % 100) == 0) report_vae(wlState.TotalSweeps, wlState.lnwlf);

        MPI_Allreduce(wlState.wlHi.data(), wlState.wlHd.data(),
                      wlState.D1BINS, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < wlState.D1BINS; ++i) {
            wlState.wlH[i] += wlState.wlHd[i];
            wlState.wllng_prior[i] += wlState.wlHd[i] * wlState.lnwlf;
        }
        std::fill(wlState.wllngi.begin(), wlState.wllngi.end(), 0.0);
        std::fill(wlState.wlHi.begin(), wlState.wlHi.end(), 0);
        tmp_flat = flatWL(WLprior);

        if ((wlState.IterSweeps % 100) == 0 && mpiState.myrank == 0) write_DOS_H();
    }

    if (fp_time) {
        std::fprintf(fp_time, "infers/ms: %f\n", 2.0 * cnt / duration);
        std::fflush(fp_time);
        std::fclose(fp_time);
    }

    FILE* fp = nullptr;
    if (mpiState.myrank == 0) fp = std::fopen("global-update-prior.dat", "w");
    for (int i = 0; i < wlState.D1BINS; ++i) {
        wlState.wllng_prior[i] *= wlState.lnwlf;
        wlState.wllng[i] += wlState.wllng_prior[i];
        if (fp) {
            std::fprintf(fp, "%g\t%18.10e\n",
                         i / wlState.invdWLD1 + wlState.WLD1min,
                         wlState.wllng_prior[i]);
        }
    }
    if (fp) std::fclose(fp);

    resetWL();
}

int WangLandau(double Ei, double Ef, SamplingMode mode) {
    const int iti = static_cast<int>((Ei * alloyState.invN - wlState.WLD1min) * wlState.invdWLD1);
    const int fti = static_cast<int>((Ef * alloyState.invN - wlState.WLD1min) * wlState.invdWLD1);

    if (fti < wlState.LOWESTE + 10 && fti >= wlState.LOWESTE) {
        const int slot = fti - wlState.LOWESTE;
        if (slot >= 0 && slot < static_cast<int>(wlState.print_list.size()) && !wlState.print_list[slot]) {
            wlState.print_list[slot] = true;
            // First visit to this low-energy slot: capture it. These are the
            // rare configurations the VAE most needs and PT sampling alone
            // rarely reaches.
            if (alloyState.snapshot_lowe) write_xyz(wlState.TotalSweeps, 1);
        }
    }

    wlState.attemptrot[iti] += 1;
    assert(iti >= 0 && iti < wlState.D1BINS);

    if (fti < 0 || fti >= wlState.D1BINS || wlState.mask[fti] == 0) {
        wlState.wlHi[iti] += 1;
        wlState.wllngi[iti] += wlState.lnwlf;
        return 0;
    }

    double lngi = 0.0, lngf = 0.0;
    if (mode == WLdos || mode == WLproduction) {
        lngi = wlState.wllngi[iti] + wlState.wllng[iti];
        lngf = wlState.wllngi[fti] + wlState.wllng[fti];
    } else {  // WLprior
        lngi = wlState.wllngi[iti] + wlState.wllng[iti] + wlState.wllng_prior[iti];
        lngf = wlState.wllngi[fti] + wlState.wllng[fti] + wlState.wllng_prior[fti];
    }

    const double R = std::exp(lngi - lngf);
    if (randd1() < R) {
        if (mode == WLprior || mode == WLdos) {
            wlState.acceptrot[iti] += 1;
        } else if (mode == WLproduction) {
            OrderParameter(fti);
        }
        wlState.wlHi[fti] += 1;
        wlState.wllngi[fti] += wlState.lnwlf;
        return 1;
    }

    if (mode == WLproduction) OrderParameter(iti);
    wlState.wlHi[iti] += 1;
    wlState.wllngi[iti] += wlState.lnwlf;
    return 0;
}
