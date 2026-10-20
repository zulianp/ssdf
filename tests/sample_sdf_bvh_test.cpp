// sample_sdf_bvh against the exact signed distance of a voxel staircase: the boundary faces of a set
// of unit lattice cells, two outward triangles each, watertight and non-convex (reflex edges), its
// faces on the sampling grid's planes and its edges on its lines, where axis-aligned rays pass
// through edges and vertices.

#include "bvh/bvh.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using G = float;
using T = float;
using I = int;

namespace {
    const int    N = 4;        // lattice cells per axis
    const double H = 0.25;     // cell size
    bool filled(const int i, const int j, const int k) {
        if (i < 0 || j < 0 || k < 0 || i >= N || j >= N || k >= N) return false;
        return i <= j && k >= 1 && k <= 2;  // a staircase in x-y, two cells thick in z
    }

    struct Face {
        int    axis;   // normal axis
        int    sign;   // outward along it
        double lo[3], hi[3];
    };
}  // namespace

int main() {
    // The lattice nodes, (N + 1)^3, and two outward triangles per boundary face.
    const int  M  = N + 1;
    auto       id = [](int i, int j, int k) { return i + M * (j + M * k); };
    std::vector<G> sx, sy, sz;
    for (int k = 0; k < M; k++) {
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                sx.push_back(G(i * H));
                sy.push_back(G(j * H));
                sz.push_back(G(k * H));
            }
        }
    }

    std::vector<I>    s0, s1, s2;
    std::vector<Face> faces;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                if (!filled(i, j, k)) continue;
                for (int axis = 0; axis < 3; axis++) {
                    for (int sign = -1; sign <= 1; sign += 2) {
                        int n[3] = {i, j, k};
                        n[axis] += sign;
                        if (filled(n[0], n[1], n[2])) continue;
                        // The face's corners, counter-clockwise seen from outside.
                        const int u = (axis + 1) % 3, v = (axis + 2) % 3;
                        int       c[4][3];
                        for (int q = 0; q < 4; q++) {
                            c[q][0] = i;
                            c[q][1] = j;
                            c[q][2] = k;
                            c[q][axis] += sign > 0 ? 1 : 0;
                        }
                        c[1][u] += 1;
                        c[2][u] += 1;
                        c[2][v] += 1;
                        c[3][v] += 1;
                        int q[4] = {0, 1, 2, 3};
                        if (sign < 0) std::swap(q[1], q[3]);
                        const int a = id(c[q[0]][0], c[q[0]][1], c[q[0]][2]), b = id(c[q[1]][0], c[q[1]][1], c[q[1]][2]),
                                  d = id(c[q[2]][0], c[q[2]][1], c[q[2]][2]), e = id(c[q[3]][0], c[q[3]][1], c[q[3]][2]);
                        s0.push_back(a), s1.push_back(b), s2.push_back(d);
                        s0.push_back(a), s1.push_back(d), s2.push_back(e);
                        Face f;
                        f.axis = axis;
                        f.sign = sign;
                        for (int dd = 0; dd < 3; dd++) {
                            f.lo[dd] = std::min(c[0][dd], c[2][dd]) * H;
                            f.hi[dd] = std::max(c[0][dd], c[2][dd]) * H;
                        }
                        faces.push_back(f);
                    }
                }
            }
        }
    }

    // The grid: spacing H/5 on [-H, (N + 1) H], so every face lies on a grid plane.
    const ptrdiff_t n  = 5 * (N + 2) + 1;
    const G         lo = G(-H), hi = G((N + 1) * H);
    std::vector<T>  out(n * n * n);
    const int rc = ssdf::sample_sdf_bvh<G, T, I>(ptrdiff_t(s0.size()), s0.data(), s1.data(), s2.data(), ptrdiff_t(sx.size()),
                                                 sx.data(), sy.data(), sz.data(), n, n, n, lo, lo, lo, hi, hi, hi, out.data());
    assert(rc == 0);

    const double h = (hi - lo) / double(n - 1);
    int          wrong_sign = 0, checked = 0;
    double       worst = 0;
    for (ptrdiff_t k = 0; k < n; k++) {
        for (ptrdiff_t j = 0; j < n; j++) {
            for (ptrdiff_t i = 0; i < n; i++) {
                const double p[3] = {lo + i * h, lo + j * h, lo + k * h};
                double       d2   = 1e300;
                for (const auto &f : faces) {
                    double s = 0;
                    for (int dd = 0; dd < 3; dd++) {
                        const double c = std::min(std::max(p[dd], f.lo[dd]), f.hi[dd]);
                        s += (p[dd] - c) * (p[dd] - c);
                    }
                    d2 = std::min(d2, s);
                }
                const double d = std::sqrt(d2);
                if (d < 1e-6) continue;  // on the surface: no sign
                const int  ci = int(std::floor(p[0] / H)), cj = int(std::floor(p[1] / H)), ck = int(std::floor(p[2] / H));
                const bool inside = filled(ci, cj, ck);
                const double got  = out[i + n * (j + n * k)];
                checked++;
                wrong_sign += (got < 0) != inside;
                worst = std::max(worst, std::abs(std::abs(got) - d));
            }
        }
    }
    std::printf("sample_sdf_bvh: %d of %d grid nodes off the surface with the wrong sign, worst distance error %g\n",
                wrong_sign,
                checked,
                worst);
    assert(wrong_sign == 0);
    assert(worst < 1e-5);
    return 0;
}
