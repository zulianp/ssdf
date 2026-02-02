#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "ssdf.hpp"

int main() {
    using G = double;
    using T = double;
    using I = int;
    std::mt19937 rng(1);
    std::uniform_real_distribution<G> dist(-1.0, 1.0);

    const int ntri = 50;
    const int npt = 20;

    std::vector<G> sx(ntri * 3), sy(ntri * 3), sz(ntri * 3);
    for (int i = 0; i < ntri * 3; ++i) {
        sx[i] = dist(rng);
        sy[i] = dist(rng);
        sz[i] = dist(rng);
    }

    std::vector<I> s0(ntri), s1(ntri), s2(ntri);
    for (int t = 0; t < ntri; ++t) {
        s0[t] = t * 3;
        s1[t] = t * 3 + 1;
        s2[t] = t * 3 + 2;
    }

    std::vector<G> x(npt), y(npt), z(npt);
    for (int i = 0; i < npt; ++i) {
        x[i] = dist(rng);
        y[i] = dist(rng);
        z[i] = dist(rng);
    }

    std::vector<T> out_ref(npt, 1e9), out_cell(npt, 1e9);

    ssdf::edf<G, T, I>(npt,
                       x.data(),
                       y.data(),
                       z.data(),
                       ntri,
                       s0.data(),
                       s1.data(),
                       s2.data(),
                       ntri * 3,
                       sx.data(),
                       sy.data(),
                       sz.data(),
                       out_ref.data());
    ssdf::edf_select<G, T, I>(npt,
                                x.data(),
                                y.data(),
                                z.data(),
                                ntri,
                                s0.data(),
                                s1.data(),
                                s2.data(),
                                ntri * 3,
                                sx.data(),
                                sy.data(),
                                sz.data(),
                                out_cell.data());

    auto brute = [&](int p) {
        G px = x[p], py = y[p], pz = z[p];
        G best = std::numeric_limits<G>::max();
        for (int t = 0; t < ntri; ++t) {
            int i0 = s0[t], i1 = s1[t], i2 = s2[t];
            G d2 = ssdf::point_triangle_dist_sq(
                px, py, pz, sx[i0], sy[i0], sz[i0], sx[i1], sy[i1], sz[i1], sx[i2], sy[i2], sz[i2]);
            best = std::min(best, d2);
        }
        return std::sqrt(best);
    };

    for (int i = 0; i < npt; ++i) {
        double brute_d = brute(i);
        {
            double d = std::abs(out_ref[i] - brute_d);
            if (d > 1e-8) {
                std::cout << "ref mismatch at " << i << " ref " << out_ref[i] << " expected " << brute_d << "\n";
                return 1;
            }
        }

        {
            double d = std::abs(out_cell[i] - brute_d);
            if (d > 1e-8) {
                std::cout << "cell mismatch at " << i << " actual " << out_ref[i] << " expected " << brute_d << "\n";
                return 1;
            }
        }
    }
    std::cout << "ref and cell matches brute" << std::endl;

    return 0;
}
