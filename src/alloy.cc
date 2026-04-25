#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include "alloy.hpp"
#include "backend/inference_backend.hpp"
#include "diagnostic.hpp"
#include "pt.hpp"
#include "rand.hpp"
#include "wanglandau.hpp"

namespace {
constexpr double kEscaleMeV = 13605.69301;
constexpr double kTscaleMeV = 0.08618;
}  // namespace

// Re-use these rather than peppering ::deepthermo_sim.constants accesses;
// the values are effectively physical constants.

void initialize() {
    alloyState.attd = alloyState.accd = 0;

    ini_coupling();
    alloyState.invN = 1.0 / (1.0 * alloyState.N * alloyState.N * alloyState.N);

    alloyState.SHIFT = alloyState.N - 1;
    alloyState.VAE_D = 2 * (alloyState.N - 1) + alloyState.SHIFT + 1;
    alloyState.PAD = static_cast<int>((std::ceil(1.0 * alloyState.VAE_D / 16) * 16 - alloyState.VAE_D) / 2);
    alloyState.VAE_D += alloyState.PAD * 2;

    const std::size_t vae_voxels =
        static_cast<std::size_t>(alloyState.VAE_D) * alloyState.VAE_D * alloyState.VAE_D * alloyState.NE;
    alloyState.inputConfig.assign(vae_voxels, 0.0f);
}

void ini_coupling() {
    const int NE = alloyState.NE;
    const int SH = alloyState.SH;
    const int n_pairs = NE * (NE - 1) / 2 + NE;

    alloyState.resize_pair_arrays();

    FILE* fop = std::fopen(alloyState.coupling_file.c_str(), "r");
    if (fop == nullptr) {
        deepthermo::die("could not open coupling file: %s",
                        alloyState.coupling_file.c_str());
    }

    std::vector<int> ti(n_pairs), tj(n_pairs);
    std::vector<double> Jr(n_pairs);
    double r = 0.0;
    double Jrold = 0.0;
    int shell = 0;
    alloyState.Nneighbors = 0;

    while (shell < SH) {
        std::fscanf(fop, "%lg", &r);
        for (int i = 0; i < n_pairs; ++i) {
            std::fscanf(fop, "%lg %d %d", &Jr[i], &ti[i], &tj[i]);
        }
        int x = 0, y = 0, z = 0;
        std::fscanf(fop, "%d %d %d\n", &x, &y, &z);

        if (alloyState.Nneighbors == 0) {
            Jrold = Jr[0];
        } else {
            if (std::fabs(r - alloyState.Dist[shell]) > 1e-6) {
                shell++;
                Jrold = Jr[0];
                if (shell == SH) break;
            }
        }
        deepthermo::Neighbor nb{x, y, z};
        alloyState.nlist.push_back(nb);
        for (int i = 0; i < n_pairs; ++i) {
            alloyState.J_at(ti[i] - 1, tj[i] - 1, shell) = Jr[i];
            alloyState.J_at(tj[i] - 1, ti[i] - 1, shell) = Jr[i];
        }
        alloyState.Dist[shell] = r;
        alloyState.NS[shell]++;
        alloyState.Nneighbors++;
    }

    std::fclose(fop);

    if (mpiState.myrank == 0) {
        for (int i = 0; i < SH; ++i) {
            std::fprintf(stderr, "%d\t%d\t%lf\t%lf\n",
                         i, alloyState.NS[i], alloyState.Dist[i], alloyState.J_at(0, 1, i));
        }
    }

    alloyState.inputPos.assign(alloyState.N_3, std::vector<int>(alloyState.Nneighbors, 0));
    for (int t = 0; t < NE; ++t) {
        alloyState.elist[t].assign(alloyState.N_3, 0);
    }
}

static void shuffle(int* array, std::size_t n) {
    for (std::size_t i = n - 1; i > 0; --i) {
        std::size_t j = std::rand() % (i + 1);
        int t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

void ini_W() {
    const int NE = alloyState.NE;
    const int SH = alloyState.SH;
    std::fill(alloyState.W.begin(), alloyState.W.end(), 0);

    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                const int idx = i * alloyState.N_2 + j * alloyState.N + k;
                const int ai = alloyState.Atom[idx];
                int cnt = 0;
                for (int shell = 0; shell < SH; ++shell) {
                    for (int ii = 0; ii < alloyState.NS[shell]; ++ii) {
                        const int aj = alloyState.Atom[alloyState.inputPos[idx][cnt++]];
                        if (ai <= aj) {
                            alloyState.W_at(ai, aj, shell)++;
                        } else {
                            alloyState.W_at(aj, ai, shell)++;
                        }
                    }
                }
                (void)NE;  // suppressed: used only for accessor indexing
            }
        }
    }
}

void ini_apos() {
    const int SH = alloyState.SH;
    std::vector<int> nn(alloyState.max_neighbors * 3);

    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                neighbor(i, j, k, nn.data());
                int cnt = 0, cntt = 0;
                const int idx = i * alloyState.N_2 + j * alloyState.N + k;
                for (int shell = 0; shell < SH; ++shell) {
                    for (int ii = 0; ii < alloyState.NS[shell]; ++ii) {
                        const int x = nn[cnt++];
                        const int y = nn[cnt++];
                        const int z = nn[cnt++];
                        alloyState.inputPos[idx][cntt++] =
                            x * alloyState.N_2 + y * alloyState.N + z;
                    }
                }
            }
        }
    }
    ini_W();
}

void ini_alloy(int state) {
    const int NE = alloyState.NE;
    const int total = alloyState.N * alloyState.N * alloyState.N;
    std::vector<int> list(total, 0);

    if (static_cast<int>(alloyState.composition.size()) != NE) {
        deepthermo::die("ini_alloy: composition has %zu entries but NE=%d",
                        alloyState.composition.size(), NE);
    }

    int cnt = 0;
    for (int t = 0; t < NE; ++t) {
        const double p = alloyState.composition[t];
        const int Ni = static_cast<int>(total * p);
        for (int cnti = 0; cnti < Ni; ++cnti) list[cnt++] = t;
        alloyState.NT[t] = 0;
    }
    while (cnt < total) {
        const int t = static_cast<int>(randd1() * NE);
        list[cnt++] = t;
    }
    if (state == 0) {
        shuffle(list.data(), total);
    }

    cnt = 0;
    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                const int t = list[cnt++];
                alloyState.Atom[i * alloyState.N_2 + j * alloyState.N + k] =
                    static_cast<short>(t);
                alloyState.NT[t]++;
            }
        }
    }

    for (int t = 0; t < NE; ++t) {
        std::fprintf(stderr, "%s:%d\n",
                     alloyState.element[t + 1].c_str(), alloyState.NT[t]);
    }
}

static inline void noffset(int i, int j, int k, int offi, int offj, int offk, int* nn, int cnt) {
    int si = i + offi;
    if (si < 0) si += alloyState.N;
    if (si >= alloyState.N) si -= alloyState.N;
    int sj = j + offj;
    if (sj < 0) sj += alloyState.N;
    if (sj >= alloyState.N) sj -= alloyState.N;
    int sk = k + offk;
    if (sk < 0) sk += alloyState.N;
    if (sk >= alloyState.N) sk -= alloyState.N;
    nn[cnt++] = si;
    nn[cnt++] = sj;
    nn[cnt++] = sk;
}

void neighbor(int i, int j, int k, int* nn) {
    int cnt = 0;
    for (int ii = 0; ii < alloyState.Nneighbors; ++ii) {
        noffset(i, j, k,
                alloyState.nlist[ii].x, alloyState.nlist[ii].y, alloyState.nlist[ii].z,
                nn, cnt);
        cnt += 3;
    }
}

static void updateWsite(int ai, int aj, int pi, int pj) {
    const int SH = alloyState.SH;
    int cnt = 0;
    for (int shell = 0; shell < SH; ++shell) {
        for (int j = 0; j < alloyState.NS[shell]; ++j) {
            const int id = alloyState.inputPos[pi][cnt++];
            if (id != pj) {
                const int a = alloyState.Atom[id];
                if (ai <= a) alloyState.W_at(ai, a, shell) -= 2;
                else         alloyState.W_at(a, ai, shell) -= 2;
                if (aj <= a) alloyState.W_at(aj, a, shell) += 2;
                else         alloyState.W_at(a, aj, shell) += 2;
            }
        }
    }
}

static void updateW(int ai, int aj, int idxi, int idxj) {
    updateWsite(ai, aj, idxi, idxj);
    updateWsite(aj, ai, idxj, idxi);
}

double Etot() {
    const int NE = alloyState.NE;
    const int SH = alloyState.SH;
    double E = 0.0;
    for (int shell = 0; shell < SH; ++shell) {
        for (int i = 0; i < NE - 1; ++i) {
            for (int j = i + 1; j < NE; ++j) {
                E += (1.0 * alloyState.W_at(i, j, shell) /
                      alloyState.NS[shell] / alloyState.N_3) *
                     alloyState.J_at(i, j, shell);
            }
        }
    }
    E += alloyState.reglin_intercept;
    return E * alloyState.N_3;
}

void BondSwap(SamplingMode mode) {
    const int i = static_cast<int>(randd1() * alloyState.N_3);
    const int n = static_cast<int>(randd1() * alloyState.Nneighbors);
    const int j = alloyState.inputPos[i][n];

    if (alloyState.Atom[i] == alloyState.Atom[j]) return;

    const double E1 = alloyState.currEtot;
    const int ai = alloyState.Atom[i];
    const int aj = alloyState.Atom[j];
    std::vector<int> Wo = alloyState.W;  // snapshot
    updateW(ai, aj, i, j);
    alloyState.Atom[i] = static_cast<short>(aj);
    alloyState.Atom[j] = static_cast<short>(ai);
    const double E2 = Etot();
    const double deltaE = E2 - E1;
    alloyState.attd++;

    auto revert = [&]() {
        alloyState.Atom[i] = static_cast<short>(ai);
        alloyState.Atom[j] = static_cast<short>(aj);
        alloyState.W = std::move(Wo);
    };

    if (mode == metropolis) {
        if (Metropolis(alloyState.currEtot, alloyState.currEtot + deltaE) == 1) {
            alloyState.currEtot += deltaE;
            alloyState.accd++;
        } else {
            revert();
        }
    } else {
        if (WangLandau(alloyState.currEtot, alloyState.currEtot + deltaE, mode) == 1) {
            alloyState.currEtot += deltaE;
            alloyState.accd++;
        } else {
            revert();
        }
    }
}

void encode(float* z) {
    const int NE = alloyState.NE;
    const int VAE_D = alloyState.VAE_D;
    const int SHIFT = alloyState.SHIFT;
    const int PAD = alloyState.PAD;

    auto voxel_index = [&](int i, int j, int k) {
        return (j - i + k + SHIFT + PAD) * VAE_D * VAE_D * NE +
               (k - j + i + SHIFT + PAD) * VAE_D * NE +
               (j + i - k + SHIFT + PAD) * NE;
    };

    std::fill(alloyState.inputConfig.begin(), alloyState.inputConfig.end(), 0.0f);
    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                const int idx = voxel_index(i, j, k);
                const int site = i * alloyState.N_2 + j * alloyState.N + k;
                const int t = alloyState.Atom[site];
                alloyState.inputConfig[idx + t] = 1.0f;
            }
        }
    }
    deepthermo_sim.backend->encode(alloyState.inputConfig.data(), VAE_D, NE, z);
}

void walk(float* npos) {
    double v1, v2, v3, a;
    v1 = 2.0 * randd1() - 1.0;
    v2 = 2.0 * randd1() - 1.0;
    a = v1 * v1 + v2 * v2;
    while (a > 1.0) {
        v1 = 2.0 * randd1() - 1.0;
        v2 = 2.0 * randd1() - 1.0;
        a = v1 * v1 + v2 * v2;
    }
    v3 = 1 - 2.0 * a;
    a = std::sqrt(1 - a);
    v1 = 2.0 * a * v1;
    v2 = 2.0 * a * v2;

    npos[0] += alloyState.Z_R * v1;
    npos[1] += alloyState.Z_R * v2;
    npos[2] += alloyState.Z_R * v3;
}

static int arg_max(const std::vector<float>& vec) {
    return static_cast<int>(std::distance(vec.begin(),
                                          std::max_element(vec.begin(), vec.end())));
}

void decode(float* z) {
    const int NE = alloyState.NE;
    const int VAE_D = alloyState.VAE_D;
    const int SHIFT = alloyState.SHIFT;
    const int PAD = alloyState.PAD;
    const std::size_t voxels =
        static_cast<std::size_t>(VAE_D) * VAE_D * VAE_D * NE;

    std::vector<int> sum(NE, 0);
    std::vector<int> sumo(NE, 0);
    std::vector<float> output(voxels, 0.0f);
    std::vector<float> vec(NE);

    auto voxel_index = [&](int i, int j, int k) {
        return (j - i + k + SHIFT + PAD) * VAE_D * VAE_D * NE +
               (k - j + i + SHIFT + PAD) * VAE_D * NE +
               (j + i - k + SHIFT + PAD) * NE;
    };

    deepthermo_sim.backend->decode(z, VAE_D, NE, output.data());

    for (int t = 0; t < NE; ++t) {
        std::fill(alloyState.elist[t].begin(), alloyState.elist[t].end(), -1);
    }
    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                const int idx = voxel_index(i, j, k);
                for (int t = 0; t < NE; ++t) vec[t] = output[idx + t];
                const int t = arg_max(vec);
                alloyState.Atom[i * alloyState.N_2 + j * alloyState.N + k] =
                    static_cast<short>(t);
                alloyState.elist[t][sum[t]++] = i * alloyState.N_2 + j * alloyState.N + k;
            }
        }
    }
    for (int t = 0; t < NE; ++t) sumo[t] = sum[t];

    // Fix concentration drift from VAE decoding (paper Algo.1 line 9: "conservation check").
    for (int t = 0; t < NE; ++t) {
        while (sum[t] > alloyState.NT[t]) {
            const int id = static_cast<int>(randd1() * sumo[t]);
            if (alloyState.elist[t][id] != -1) {
                int tid = -1;
                for (int tt = 0; tt < NE; ++tt) {
                    if (tt != t && sum[tt] < alloyState.NT[tt]) {
                        tid = tt;
                        break;
                    }
                }
                alloyState.Atom[alloyState.elist[t][id]] = static_cast<short>(tid);
                alloyState.elist[t][id] = -1;
                sum[t]--;
                sum[tid]++;
            }
        }
    }
    for (int t = 0; t < NE; ++t) assert(sum[t] == alloyState.NT[t]);
    ini_W();
}

void vae_update(SamplingMode mode) {
    float z[3] = {0, 0, 0};
    alloyState.Atomo = alloyState.Atom;
    std::vector<int> Wo = alloyState.W;

    encode(z);
    walk(z);
    decode(z);

    const double E1 = alloyState.currEtot;
    const double E2 = Etot();
    const double deltaE = E2 - E1;
    alloyState.attd++;

    auto revert = [&]() {
        alloyState.Atom = alloyState.Atomo;
        alloyState.W = std::move(Wo);
    };

    if (mode == metropolis) {
        if (Metropolis(alloyState.currEtot, alloyState.currEtot + deltaE) == 1) {
            alloyState.currEtot += deltaE;
            alloyState.accd++;
        } else {
            revert();
        }
    } else {
        if (WangLandau(alloyState.currEtot, alloyState.currEtot + deltaE, mode) == 1) {
            alloyState.currEtot += deltaE;
            alloyState.accd++;
        } else {
            revert();
        }
    }
}

void write_pos() {
    FILE* ofp = std::fopen("compos.dat", "w");
    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                std::fprintf(ofp, "%d\t%d\t%d\t%d\t \n", i + k, i + j, j + k,
                             alloyState.Atom[i * alloyState.N_2 + j * alloyState.N + k]);
            }
        }
    }
    std::fclose(ofp);
}

void write_xyz(int frame) {
    const int NE = alloyState.NE;
    const int SH = alloyState.SH;
    char s[512];
    std::snprintf(s, sizeof(s), "snap_%d_%d.xyz", frame, mpiState.myrank);
    FILE* ofp = std::fopen(s, "ab+");
    std::fprintf(ofp, "%d\n", alloyState.N_3);
    std::fprintf(ofp, "Eng = %.8lg  MC_step = %d\n", alloyState.currEtot, frame);

    for (int shell = 0; shell < SH; ++shell) {
        for (int i = 0; i < NE - 1; ++i) {
            for (int j = i + 1; j < NE; ++j) {
                std::fprintf(ofp, "%.8lg\t",
                             1.0 * alloyState.W_at(i, j, shell) /
                                 alloyState.NS[shell] / alloyState.N_3);
            }
        }
    }
    std::fprintf(ofp, "\n");

    for (int i = 0; i < alloyState.N; ++i) {
        for (int j = 0; j < alloyState.N; ++j) {
            for (int k = 0; k < alloyState.N; ++k) {
                std::fprintf(ofp, "%s %d %d %d \n",
                             alloyState.element[alloyState.Atom[i * alloyState.N_2 + j * alloyState.N + k] + 1].c_str(),
                             j - i + k, k - j + i, j + i - k);
            }
        }
    }
    std::fclose(ofp);
}

void thermoqs() {
    FILE* therm_op = std::fopen("therm.dat", "w");
    if (therm_op == nullptr) {
        deepthermo::die("thermoqs(): could not open therm.dat");
    }
    double Nfree = 0.0;

    for (double T = alloyState.TTi; T <= alloyState.TTf + alloyState.dTT;
         T += alloyState.dTT) {
        double U = 0.0, Z = 0.0, C = 0.0, F = 0.0, S = 0.0;
        double Bw = 0.0;
        double M = 0.0, M2 = 0.0, X = 0.0;
        double lambda = -1.0e300;
        double centerE = 0.0;

        for (int i = 0; i < wlState.D1BINS; ++i) {
            centerE = ((i / wlState.invdWLD1 + wlState.WLD1min) +
                       0.5 * (wlState.WLD1max - wlState.WLD1min) / (1.0 * wlState.D1BINS)) /
                      alloyState.invN;
            if (lambda < (wlState.wllng[i] - 1.0 * centerE * kEscaleMeV / (T * kTscaleMeV))) {
                lambda = wlState.wllng[i] - 1.0 * centerE * kEscaleMeV / (T * kTscaleMeV);
            }
        }

        for (int i = 0; i < wlState.D1BINS; ++i) {
            centerE = ((i / wlState.invdWLD1 + wlState.WLD1min) +
                       0.5 * (wlState.WLD1max - wlState.WLD1min) / (1.0 * wlState.D1BINS)) /
                      alloyState.invN;
            Bw = std::exp(wlState.wllng[i] - centerE * kEscaleMeV / (T * kTscaleMeV) - lambda);
            Z += Bw;
            U += centerE * Bw;
            C += centerE * centerE * Bw;
            if (wlState.wlH[i] > 0) {
                M += wlState.wllngi[i] / wlState.wlH[i] * Bw;
                M2 += wlState.wllngd[i] / wlState.wlH[i] * Bw;
            }
        }

        U = U / Z;
        if (T == alloyState.TTi) {
            Nfree = -(T * kTscaleMeV) / kEscaleMeV * (lambda + std::log(Z)) - U;
        }
        C = ((C / Z) - (U * U)) * kEscaleMeV * kEscaleMeV /
            (T * T * kTscaleMeV * kTscaleMeV);
        F = -(T * kTscaleMeV) / kEscaleMeV * (lambda + std::log(Z)) - Nfree;
        S = (U - F) * kEscaleMeV / (T * kTscaleMeV);
        M /= Z;
        M2 /= Z;
        X = (M2 - M * M) * alloyState.N_3 / T / kTscaleMeV;

        std::fprintf(therm_op, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                     T, U * alloyState.invN, C * alloyState.invN, F * alloyState.invN,
                     S * alloyState.invN, M, X);
    }
    std::fclose(therm_op);
}

double L1() {
    float z[3];
    encode(z);
    double M = 0.0;
    for (int i = 0; i < 3; ++i) M += std::fabs(z[i]);
    return M;
}

void OrderParameter(int idx) {
    const double M = L1();
    assert(idx >= 0 && idx < wlState.D1BINS);
    alloyState.op[idx] += M;
    alloyState.op2[idx] += M * M;
}
