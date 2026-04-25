#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>
#include <vector>

#include "mpi.h"

#include "alloy.hpp"
#include "diagnostic.hpp"
#include "pt.hpp"
#include "rand.hpp"

namespace {

constexpr double kEscaleMeV = 13605.69301;
constexpr double kTscaleMeV = 0.08618;

// Per-process Mersenne Twister, seeded once in ini_sys(). Replaces the
// header-defined random.h that used to live alongside this TU.
std::mt19937& mt_engine() {
    static std::mt19937 engine{std::random_device{}()};
    return engine;
}

void mt_seed(unsigned long s) { mt_engine().seed(s); }
double mt_uniform()           { return std::uniform_real_distribution<double>{0.0, 1.0}(mt_engine()); }

}  // namespace

void ini_T(double Ti, double Tf, int nT) {
    if (Ti * Tf <= 1e-5) {
        deepthermo::die("ini_T: invalid temperature range (Ti=%g, Tf=%g)", Ti, Tf);
    }
    ptState.T.assign(nT, 0.0);
    ptState.T[0] = std::log(Ti);
    ptState.T[nT - 1] = std::log(Tf);
    for (int i = 1; i < nT - 1; ++i) {
        ptState.T[i] = ptState.T[0] + (ptState.T[nT - 1] - ptState.T[0]) * i / (nT - 1);
    }
    for (int i = 0; i < nT; ++i) ptState.T[i] = std::exp(ptState.T[i]);
    ptState.pT = ptState.T[mpiState.myrank];
}

void ini_sys() {
    alloyState.N_2 = alloyState.N * alloyState.N;
    alloyState.N_3 = alloyState.N * alloyState.N * alloyState.N;

    alloyState.Atom.assign(alloyState.N_3, 0);
    alloyState.Atomo.assign(alloyState.N_3, 0);
    ptState.att = ptState.acc = 0;
    mt_seed(static_cast<unsigned long>(mpiState.myrank * std::rand()));

    initialize();
}

namespace {

void send_config_swap(int peer, short* buf_ui, short* buf_uj, std::size_t mem_u) {
    MPI_Status istatus;
    std::memcpy(buf_ui, alloyState.Atom.data(), mem_u);
    MPI_Sendrecv(buf_ui, alloyState.N_3, MPI_SHORT, peer, 123,
                 buf_uj, alloyState.N_3, MPI_SHORT, peer, 123,
                 MPI_COMM_WORLD, &istatus);
    std::memcpy(alloyState.Atom.data(), buf_uj, mem_u);
}

}  // namespace

void swap(bool even) {
    constexpr int itag = 123;
    MPI_Status istatus;
    const std::size_t mem_u = alloyState.N_3 * sizeof(short);
    std::vector<short> buf_ui(alloyState.N_3);
    std::vector<short> buf_uj(alloyState.N_3);

    const int peer_up = mpiState.myrank + 1;
    const int peer_dn = mpiState.myrank - 1;

    auto propose_up = [&]() {
        double Ei = alloyState.currEtot;
        double Ej = 0.0;
        MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, peer_up, itag,
                     &Ej, 1, MPI_DOUBLE, peer_up, itag, MPI_COMM_WORLD, &istatus);
        const double delta = (Ei - Ej) * kEscaleMeV *
                             (1 / (kTscaleMeV * ptState.T[peer_up]) -
                              1 / (kTscaleMeV * ptState.T[mpiState.myrank]));
        bool flag = (delta <= 0 || mt_uniform() < std::exp(-delta));
        MPI_Send(&flag, 1, MPI_LOGICAL, peer_up, itag, MPI_COMM_WORLD);
        ptState.att++;
        if (flag) {
            send_config_swap(peer_up, buf_ui.data(), buf_uj.data(), mem_u);
            alloyState.currEtot = Ej;
            ptState.acc++;
        }
    };

    auto accept_from_dn = [&]() {
        double Ei = alloyState.currEtot;
        double Ej = 0.0;
        MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, peer_dn, itag,
                     &Ej, 1, MPI_DOUBLE, peer_dn, itag, MPI_COMM_WORLD, &istatus);
        bool flag = false;
        MPI_Recv(&flag, 1, MPI_LOGICAL, peer_dn, itag, MPI_COMM_WORLD, &istatus);
        if (flag) {
            send_config_swap(peer_dn, buf_ui.data(), buf_uj.data(), mem_u);
            alloyState.currEtot = Ej;
        }
    };

    if (even) {
        if (mpiState.myrank % 2 == 0) {
            propose_up();
        } else {
            accept_from_dn();
        }
    } else {
        if (mpiState.myrank % 2 != 0 && mpiState.myrank != mpiState.nprocs - 1) {
            propose_up();
        } else if (mpiState.myrank != 0 && mpiState.myrank != mpiState.nprocs - 1) {
            accept_from_dn();
        }
    }

    ini_W();
}

int Metropolis(double Ei, double Ef) {
    if (Ef < Ei) return 1;
    const double R = std::exp(-(Ef - Ei) * kEscaleMeV / (kTscaleMeV * ptState.pT));
    return (randd1() < R) ? 1 : 0;
}

void mchybrid(SamplingMode mode) {
    const int total = alloyState.N * alloyState.N * alloyState.N;
    for (int cnt = 0; cnt < total; ++cnt) BondSwap(mode);
}

void parallel_tempering(int nT, double DROPI, double SAMPS, double SEP, int irun,
                        SamplingMode mode) {
    std::vector<double> Ea(nT, 0.0);
    std::vector<double> Ma(nT, 0.0);
    std::vector<double> ca(nT, 0.0);
    std::vector<double> xa(nT, 0.0);
    std::vector<double> bca(nT, 0.0);
    std::vector<int> attempts(nT, 0);
    std::vector<int> accepts(nT, 0);

    double avgE = 0.0, avgE2 = 0.0, avgM = 0.0, avgM2 = 0.0, avgM4 = 0.0;
    bool flag = true, stop = false;
    std::time_t t1 = 0, t2 = 0;

    if (mpiState.myrank == 0) t1 = std::time(nullptr);

    if (ptState.Restart == 0) {
        for (int mcs = 0; mcs < DROPI; ++mcs) {
            mchybrid(mode);
            if (mcs % 2 == 0) {
                swap(flag);
                flag = !flag;
            }
        }
    }

    int mcs = 0;
    if (ptState.Restart == 1) {
        char s[512];
        std::snprintf(s, sizeof(s), "mc%d.input", mpiState.myrank);
        FILE* ofp1 = std::fopen(s, "rb");
        if (ofp1 == nullptr) deepthermo::die("mc.input not found: %s", s);
        std::fread(&mcs, sizeof(int), 1, ofp1);
        std::fread(&avgE, sizeof(double), 1, ofp1);
        std::fread(&avgE2, sizeof(double), 1, ofp1);
        std::fread(&avgM, sizeof(double), 1, ofp1);
        std::fread(&avgM2, sizeof(double), 1, ofp1);
        std::fread(&avgM4, sizeof(double), 1, ofp1);
        std::fclose(ofp1);
    }

    while (mcs < SAMPS) {
        for (int i = 0; i < SEP; ++i) {
            mchybrid(mode);
            if (static_cast<int>(mcs * SEP + i) % 2 == 0) {
                swap(flag);
                flag = !flag;
            }
        }

        avgE += alloyState.currEtot;
        avgE2 += alloyState.currEtot * alloyState.currEtot;
        const double M = L1();
        avgM += M;
        avgM2 += M * M;
        avgM4 += M * M;

        if (mcs % CHPT_STEPS == 0) {
            std::printf("mpiState.myrank =%d, currEtot = %g, Eng = %g\n",
                        mpiState.myrank, alloyState.currEtot, Etot());
            if (mpiState.myrank == 0) {
                t2 = std::time(nullptr);
                if (1.0 * (t2 - t1) / 3600 > 0.9 * TIMER) stop = true;
            }
            MPI_Bcast(&stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
            if (stop) {
                write_state();
                char s[512];
                std::snprintf(s, sizeof(s), "mc%d.input", mpiState.myrank);
                FILE* ofp1 = std::fopen(s, "wb");
                mcs++;
                std::fwrite(&mcs, sizeof(int), 1, ofp1);
                std::fwrite(&avgE, sizeof(double), 1, ofp1);
                std::fwrite(&avgE2, sizeof(double), 1, ofp1);
                std::fwrite(&avgM, sizeof(double), 1, ofp1);
                std::fwrite(&avgM2, sizeof(double), 1, ofp1);
                std::fwrite(&avgM4, sizeof(double), 1, ofp1);
                std::fclose(ofp1);
                std::printf("avgE2-avgE*avgE %g\n", avgE2 / mcs - avgE * avgE / mcs / mcs);
                MPI_Finalize();
                std::exit(123);
            }
        }
        mcs++;
    }

    if (mpiState.myrank == 0) write_pos();

    avgE /= 1.0 * SAMPS;
    avgE2 /= 1.0 * SAMPS;
    avgM /= 1.0 * SAMPS;
    avgM2 /= 1.0 * SAMPS;
    avgM4 /= 1.0 * SAMPS;
    const double c = (avgE2 - avgE * avgE) * kEscaleMeV * kEscaleMeV /
                     (kTscaleMeV * kTscaleMeV * ptState.pT * ptState.pT);
    const double x = (avgM2 - avgM * avgM) * alloyState.N_3 / ptState.pT / kTscaleMeV;
    const double bc = 1 - avgM4 / avgM2 / avgM2 / 3;

    MPI_Gather(&avgE, 1, MPI_DOUBLE, Ea.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&avgM, 1, MPI_DOUBLE, Ma.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&c, 1, MPI_DOUBLE, ca.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&x, 1, MPI_DOUBLE, xa.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&bc, 1, MPI_DOUBLE, bca.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&ptState.att, 1, MPI_INT, attempts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&ptState.acc, 1, MPI_INT, accepts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    FILE* ofp2 = nullptr;
    if (mpiState.myrank == 0) {
        char s[512];
        std::snprintf(s, sizeof(s), "stat%d.dat", irun);
        FILE* ofp = std::fopen(s, "w");
        for (int n = 0; n < nT; ++n) {
            std::fprintf(ofp, "%g\t%g\t%g\t",
                         ptState.T[n], Ea[n] * alloyState.invN, ca[n] * alloyState.invN);
            std::fprintf(ofp, "%g\t%g\t%g\t", Ma[n], xa[n], bca[n]);
            std::fprintf(ofp, "\n");
        }
        std::fclose(ofp);

        std::snprintf(s, sizeof(s), "misc%d.dat", irun);
        ofp2 = std::fopen(s, "w");
        std::fprintf(ofp2, "swap prob:\n");
        for (int n = 0; n < nT - 1; ++n) {
            const double ratio = attempts[n] > 0 ? 1.0 * accepts[n] / attempts[n] : 0.0;
            std::fprintf(ofp2, "%g\t%g\n", ptState.T[n], ratio);
        }
    }

    MPI_Gather(&alloyState.attd, 1, MPI_INT, attempts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&alloyState.accd, 1, MPI_INT, accepts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpiState.myrank == 0 && ofp2) {
        std::fprintf(ofp2, "Rot prob:\n");
        for (int n = 0; n < nT; ++n) {
            const double ratio = attempts[n] > 0 ? 1.0 * accepts[n] / attempts[n] : 0.0;
            std::fprintf(ofp2, "%g\t%g\n", ptState.T[n], ratio);
        }
        std::fflush(ofp2);
        std::fclose(ofp2);
    }
}

void write_state() {
    char s[512];
    std::snprintf(s, sizeof(s), "state%d.input", mpiState.myrank);
    FILE* fp = std::fopen(s, "wb");
    std::fwrite(alloyState.Atom.data(), sizeof(short), alloyState.N_3, fp);
    std::fwrite(alloyState.NT.data(), sizeof(int), alloyState.NE, fp);
    std::fwrite(&alloyState.currEtot, sizeof(double), 1, fp);
    std::fclose(fp);
}

void read_state() {
    char s[512];
    std::snprintf(s, sizeof(s), "state%d.input", mpiState.myrank);
    FILE* fp = std::fopen(s, "rb");
    if (fp == nullptr) deepthermo::die("state file not found: %s", s);
    double eng = 0.0;
    std::fread(alloyState.Atom.data(), sizeof(short), alloyState.N_3, fp);
    std::fread(alloyState.NT.data(), sizeof(int), alloyState.NE, fp);
    std::fread(&eng, sizeof(double), 1, fp);
    ini_apos();
    alloyState.currEtot = Etot();
    if (std::fabs(eng - alloyState.currEtot) > 1e-5) {
        deepthermo::die("state file corrupted: stored E=%g, recomputed E=%g",
                        eng, alloyState.currEtot);
    }
    std::fclose(fp);
}
