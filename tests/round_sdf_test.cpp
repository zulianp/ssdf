// rebuild_sdf and round_sdf against exact signed distances:
//  - a tilted plane, whose level set marching tetrahedra reproduce exactly: the rebuild returns the
//    field (to float rounding) where the closest point lies inside the grid;
//  - a box rounded by r: opening the box (f > 0 inside it) and closing its complement (f > 0
//    outside) both give the rounded box, the Minkowski sum of the box shrunk by r and a ball of
//    radius r, whose distance is known; the error must fall with the grid spacing.

#include "bvh/sdf_rebuild.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using G = float;
using T = float;

namespace {

    // The signed distance to the box |x_d| <= b_d, rounded by r (r = 0: the box), positive outside.
    double round_box(const double *p, const double *b, const double r) {
        double q[3], out = 0, in = -1e300;
        for (int d = 0; d < 3; ++d) {
            q[d] = std::fabs(p[d]) - b[d] + r;
            out += std::max(q[d], 0.0) * std::max(q[d], 0.0);
            in = std::max(in, q[d]);
        }
        return std::sqrt(out) + std::min(in, 0.0) - r;
    }

    ssdf::GridGeometry<G> cube(const ptrdiff_t n) {
        ssdf::GridGeometry<G> g;
        for (int d = 0; d < 3; ++d) {
            g.n[d] = n;
            g.o[d] = G(-0.5);
            g.h[d] = G(1.0 / double(n - 1));
        }
        return g;
    }

    void point(const ssdf::GridGeometry<G> &g, const ptrdiff_t i, const ptrdiff_t j, const ptrdiff_t k, double *p) {
        p[0] = double(g.o[0]) + double(i) * double(g.h[0]);
        p[1] = double(g.o[1]) + double(j) * double(g.h[1]);
        p[2] = double(g.o[2]) + double(k) * double(g.h[2]);
    }

    // The largest error of a rounding of the box sampled at 1/(n - 1) and refined by `refine` first,
    // written back on the 1/(n - 1) samples.
    double refined_box_error(const ptrdiff_t n, const double sign, const ptrdiff_t refine) {
        const double b[3] = {0.3, 0.25, 0.2}, r = 0.08;
        const auto   g    = cube(n);
        std::vector<T> f(g.size());
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    f[i + n * (j + n * k)] = T(-sign * round_box(p, b, 0));
                }
            }
        }
        const auto     fg = g.refined(refine);
        std::vector<T> fine(fg.size()), out(g.size());
        ssdf::refine_field<G, T>(g, f.data(), refine, fine.data());
        const int rc = ssdf::round_sdf<G, T>(fg, fine.data(), T(r), T(0.3), true, true, refine, out.data(), nullptr);
        assert(rc == 0);
        (void)rc;
        double err = 0;
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    const double exact = -sign * round_box(p, b, r);
                    if (std::fabs(exact) > r) continue;
                    err = std::max(err, std::fabs(double(out[i + n * (j + n * k)]) - exact));
                }
            }
        }
        return err;
    }

    // The largest |out - sign * round_box| over the samples within r of the rounded box's surface;
    // and the largest change of out by one more rebuild (zero for an exact distance to its zero
    // level set, up to the triangles' approximation of it).
    double rounded_box_error(const ptrdiff_t n, const double sign, const ptrdiff_t stride, double *self = nullptr) {
        const double b[3] = {0.3, 0.25, 0.2}, r = 0.08;
        const auto   g    = cube(n);
        std::vector<T> f(g.size());
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    f[i + n * (j + n * k)] = T(-sign * round_box(p, b, 0));
                }
            }
        }
        const auto     og = g.strided(stride);
        std::vector<T> out(og.size());
        const int      rc = ssdf::round_sdf<G, T>(g, f.data(), T(r), T(0.3), true, true, stride, out.data(), nullptr);
        assert(rc == 0);
        (void)rc;
        std::vector<T> again(og.size());
        if (self) {
            const int rc2 = ssdf::rebuild_sdf<G, T>(og, out.data(), T(0), T(0.3), 1, again.data());
            assert(rc2 == 0);
            (void)rc2;
            *self = 0;
        }
        double err = 0;
        for (ptrdiff_t k = 0; k < og.n[2]; ++k) {
            for (ptrdiff_t j = 0; j < og.n[1]; ++j) {
                for (ptrdiff_t i = 0; i < og.n[0]; ++i) {
                    double p[3];
                    point(og, i, j, k, p);
                    const double    exact = -sign * round_box(p, b, r);
                    const ptrdiff_t q     = i + og.n[0] * (j + og.n[1] * k);
                    if (std::fabs(exact) > r) continue;
                    err = std::max(err, std::fabs(double(out[q]) - exact));
                    if (self) *self = std::max(*self, std::fabs(double(again[q]) - double(out[q])));
                }
            }
        }
        return err;
    }

}  // namespace

int main() {
    // A tilted plane: n . x = c.
    {
        const ptrdiff_t n = 33;
        const auto      g = cube(n);
        const double    nrm[3] = {0.48, 0.6, 0.64}, c = 0.05, band = 0.2;  // |nrm| = 1
        std::vector<T>  f(g.size()), out(g.size());
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    f[i + n * (j + n * k)] = T(nrm[0] * p[0] + nrm[1] * p[1] + nrm[2] * p[2] - c);
                }
            }
        }
        const int rc = ssdf::rebuild_sdf<G, T>(g, f.data(), T(0), T(band), 1, out.data());
        assert(rc == 0);
        (void)rc;
        double err = 0;
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    const double s = nrm[0] * p[0] + nrm[1] * p[1] + nrm[2] * p[2] - c;
                    // The closest point on the plane, inside the grid by a cell.
                    bool inside = true;
                    for (int d = 0; d < 3; ++d) inside &= std::fabs(p[d] - s * nrm[d]) < 0.5 - double(g.h[d]);
                    if (!inside || std::fabs(s) > band) continue;
                    err = std::max(err, std::fabs(double(out[i + n * (j + n * k)]) - s));
                }
            }
        }
        printf("plane: max error %.3g\n", err);
        assert(err < 1e-6);
    }

    // refine_field reproduces a cubic away from the grid's ends (where the end samples repeat).
    {
        const ptrdiff_t n = 9, m = 3;
        const auto      g = cube(n), fg = g.refined(m);
        auto            cubic = [](const double *p) {
            return p[0] * p[0] * p[0] - 2 * p[1] * p[1] + p[0] * p[1] * p[2] + 0.5 * p[2] - 0.1;
        };
        std::vector<T> f(g.size()), fine(fg.size());
        for (ptrdiff_t k = 0; k < n; ++k) {
            for (ptrdiff_t j = 0; j < n; ++j) {
                for (ptrdiff_t i = 0; i < n; ++i) {
                    double p[3];
                    point(g, i, j, k, p);
                    f[i + n * (j + n * k)] = T(cubic(p));
                }
            }
        }
        ssdf::refine_field<G, T>(g, f.data(), m, fine.data());
        double err = 0;
        for (ptrdiff_t k = m; k + m < fg.n[2]; ++k) {
            for (ptrdiff_t j = m; j + m < fg.n[1]; ++j) {
                for (ptrdiff_t i = m; i + m < fg.n[0]; ++i) {
                    double p[3];
                    point(fg, i, j, k, p);
                    err = std::max(err, std::fabs(double(fine[i + fg.n[0] * (j + fg.n[1] * k)]) - cubic(p)));
                }
            }
        }
        printf("refine_field: cubic reproduced to %.3g\n", err);
        assert(err < 1e-6);
    }

    // The rounded box, by opening the box and by closing its complement, at three spacings; and the
    // finest but one written on the coarsest's samples.
    for (const double sign : {1.0, -1.0}) {
        double       s1, s2, s3;
        const double e1 = rounded_box_error(41, sign, 1, &s1);
        const double e2 = rounded_box_error(81, sign, 1, &s2);
        const double e3 = rounded_box_error(161, sign, 1, &s3);
        const double es = rounded_box_error(81, sign, 2);
        const double er = refined_box_error(41, sign, 2);
        printf("rounded box (%s), h = 1/40, 1/80, 1/160: max error %.3g, %.3g, %.3g (rates %.2f, %.2f); "
               "change by a rebuild %.3g, %.3g, %.3g (rates %.2f, %.2f); 1/80 written at 1/40 %.3g; 1/40 refined by 2 %.3g\n",
               sign > 0 ? "opening the box" : "closing its complement",
               e1,
               e2,
               e3,
               std::log2(e1 / e2),
               std::log2(e2 / e3),
               s1,
               s2,
               s3,
               std::log2(s1 / s2),
               std::log2(s2 / s3),
               es,
               er);
        // The rounding converges (first order where the eroded or dilated box's edges were cut by
        // the cells, faster elsewhere); the result is a distance to its own zero level set at
        // second order; the strided output carries the finer grid's accuracy.
        assert(e1 < 1.0 / 40 && e2 < 0.6 / 80 && e3 < 0.25 / 160);
        assert(e2 < 0.6 * e1 && e3 < 0.6 * e2);
        assert(s1 < 0.1 / 40 && s2 < 0.35 * s1 && s3 < 0.35 * s2);
        assert(es <= 1.01 * e2);
        // Refined by 2, the 1/40 samples get most of the 1/80 grid's accuracy (the box's edges come
        // through the cubic smoothed, then the rounding replaces them).
        assert(er < 0.6 * e1);
    }
    return 0;
}
