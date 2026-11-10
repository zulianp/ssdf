#include "sdf_rebuild.hpp"

#include "cuBQL/builder/cpu.h"
#include "cuBQL/queries/triangleData/closestPointOnAnyTriangle.h"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>

namespace ssdf {

    namespace {

        // The six tetrahedra of a cell about its diagonal from corner 0 to corner 7; corner c is at
        // (c & 1, (c >> 1) & 1, (c >> 2) & 1) in cell units.
        constexpr int kTets[6][4] = {{0, 1, 3, 7}, {0, 1, 5, 7}, {0, 2, 3, 7}, {0, 2, 6, 7}, {0, 4, 5, 7}, {0, 4, 6, 7}};

        struct Cell {
            ptrdiff_t g[8];     // grid index
            double    v[8];     // value less the level
            double    p[8][3];  // position
        };

        // Cell (i, j, k), its corners from (i, j, k) to (i + 1, j + 1, k + 1). False when all its
        // corners lie on one side of the level set (no triangle).
        template <typename G, typename T>
        bool load_cell(const GridGeometry<G> &grid,
                       const T *const         f,
                       const T                level,
                       const ptrdiff_t        i,
                       const ptrdiff_t        j,
                       const ptrdiff_t        k,
                       Cell                  &cell) {
            const ptrdiff_t sy = grid.n[0], sz = grid.n[0] * grid.n[1];
            int             above = 0;
            for (int c = 0; c < 8; ++c) {
                const ptrdiff_t ci = i + (c & 1), cj = j + ((c >> 1) & 1), ck = k + ((c >> 2) & 1);
                cell.g[c]          = ci + cj * sy + ck * sz;
                cell.v[c]          = double(f[cell.g[c]]) - double(level);
                above += cell.v[c] > 0;
            }
            if (above == 0 || above == 8) return false;
            for (int c = 0; c < 8; ++c) {
                const ptrdiff_t ijk[3] = {i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1)};
                for (int d = 0; d < 3; ++d) {
                    cell.p[c][d] = double(grid.o[d]) + double(ijk[d]) * double(grid.h[d]);
                }
            }
            return true;
        }

        // The crossing on the edge between corners a and b (one above the level, one not), from the
        // corner of lower grid index, so that every cell sharing the edge finds the same point.
        void crossing(const Cell &cell, int a, int b, double *const q) {
            if (cell.g[b] < cell.g[a]) {
                const int s = a;
                a           = b;
                b           = s;
            }
            const double t = cell.v[a] / (cell.v[a] - cell.v[b]);
            for (int d = 0; d < 3; ++d) {
                q[d] = cell.p[a][d] + t * (cell.p[b][d] - cell.p[a][d]);
            }
        }

        // The triangles (0, 1 or 2) of tetrahedron t of the cell, their vertices into tri, without
        // those with two vertices equal in G: a sample on the level set makes crossings coincide, and
        // a zero-length edge turns cuBQL's point-to-edge test into NaN, which then wins every later
        // comparison of the query (the closest triangle is lost).
        template <typename G>
        int tet_triangles(const Cell &cell, const int *const t, double (*const tri)[3][3]) {
            int up[4], down[4], nu = 0, nd = 0;
            for (int a = 0; a < 4; ++a) {
                if (cell.v[t[a]] > 0) {
                    up[nu++] = t[a];
                } else {
                    down[nd++] = t[a];
                }
            }
            if (nu == 0 || nd == 0) return 0;
            int m = 0;
            if (nu == 1 || nd == 1) {
                // One corner apart: a triangle on its three edges.
                const int  lone   = nu == 1 ? up[0] : down[0];
                const int *others = nu == 1 ? down : up;
                for (int a = 0; a < 3; ++a) crossing(cell, lone, others[a], tri[0][a]);
                m = 1;
            } else {
                // Two and two: the quadrilateral on the four edges between them, as two triangles.
                double q[4][3];
                crossing(cell, up[0], down[0], q[0]);
                crossing(cell, up[0], down[1], q[1]);
                crossing(cell, up[1], down[1], q[2]);
                crossing(cell, up[1], down[0], q[3]);
                const int quad[2][3] = {{0, 1, 2}, {0, 2, 3}};
                for (int s = 0; s < 2; ++s) {
                    for (int a = 0; a < 3; ++a) {
                        for (int d = 0; d < 3; ++d) tri[s][a][d] = q[quad[s][a]][d];
                    }
                }
                m = 2;
            }
            int kept = 0;
            for (int s = 0; s < m; ++s) {
                bool degenerate = false;
                for (int a = 0; a < 3; ++a) {
                    const int b = (a + 1) % 3;
                    degenerate |= G(tri[s][a][0]) == G(tri[s][b][0]) && G(tri[s][a][1]) == G(tri[s][b][1]) &&
                                  G(tri[s][a][2]) == G(tri[s][b][2]);
                }
                if (degenerate) continue;
                if (kept != s) {
                    for (int a = 0; a < 3; ++a) {
                        for (int d = 0; d < 3; ++d) tri[kept][a][d] = tri[s][a][d];
                    }
                }
                ++kept;
            }
            return kept;
        }

    }  // namespace

    template <typename G>
    GridGeometry<G> GridGeometry<G>::strided(const ptrdiff_t stride) const {
        GridGeometry<G> s = *this;
        for (int d = 0; d < 3; ++d) {
            s.n[d] = (n[d] - 1) / stride + 1;
            s.h[d] = h[d] * G(stride);
        }
        return s;
    }

    template <typename G>
    GridGeometry<G> GridGeometry<G>::refined(const ptrdiff_t factor) const {
        GridGeometry<G> r = *this;
        for (int d = 0; d < 3; ++d) {
            r.n[d] = (n[d] - 1) * factor + 1;
            r.h[d] = G(double(h[d]) / double(factor));
        }
        return r;
    }

    template <typename G, typename T>
    void refine_field(const GridGeometry<G> &grid, const T *const f, const ptrdiff_t factor, T *const fine) {
        const GridGeometry<G> rg = grid.refined(factor);
        // The weights of the cubic through samples -1, 0, 1, 2 at t = m / factor in [0, 1).
        std::vector<double> w(4 * factor);
        for (ptrdiff_t m = 0; m < factor; ++m) {
            const double t   = double(m) / double(factor);
            w[4 * m + 0]     = -t * (t - 1) * (t - 2) / 6;
            w[4 * m + 1]     = (t + 1) * (t - 1) * (t - 2) / 2;
            w[4 * m + 2]     = -(t + 1) * t * (t - 2) / 2;
            w[4 * m + 3]     = (t + 1) * t * (t - 1) / 6;
        }
        // Axis by axis: x on the coarse grid into the fine one's x extent, then y, then z, in place
        // in `fine` (each pass reads the samples the previous one placed at multiples of factor).
        const ptrdiff_t fx = rg.n[0], fxy = rg.n[0] * rg.n[1];
#pragma omp parallel for collapse(2)
        for (ptrdiff_t k = 0; k < grid.n[2]; ++k) {
            for (ptrdiff_t j = 0; j < grid.n[1]; ++j) {
                for (ptrdiff_t i = 0; i < grid.n[0]; ++i) {
                    fine[i * factor + fx * (j * factor) + fxy * (k * factor)] = f[i + grid.n[0] * (j + grid.n[1] * k)];
                }
            }
        }
        for (int axis = 0; axis < 3; ++axis) {
            const ptrdiff_t stride = axis == 0 ? 1 : (axis == 1 ? fx : fxy);
            const ptrdiff_t n      = grid.n[axis];
            // The lines along `axis`: every fine index on the other axes up to their current extent,
            // at multiples of factor on the axes not yet refined.
            const int       a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
            const ptrdiff_t s1 = a1 == 0 ? 1 : (a1 == 1 ? fx : fxy), s2 = a2 == 0 ? 1 : (a2 == 1 ? fx : fxy);
            const ptrdiff_t c1 = a1 < axis ? 1 : factor, c2 = a2 < axis ? 1 : factor;
            const ptrdiff_t n1 = a1 < axis ? rg.n[a1] : grid.n[a1], n2 = a2 < axis ? rg.n[a2] : grid.n[a2];
#pragma omp parallel for collapse(2)
            for (ptrdiff_t q2 = 0; q2 < n2; ++q2) {
                for (ptrdiff_t q1 = 0; q1 < n1; ++q1) {
                    T *const line = fine + q1 * c1 * s1 + q2 * c2 * s2;
                    for (ptrdiff_t c = 0; c + 1 < n; ++c) {
                        const T v[4] = {line[std::max<ptrdiff_t>(c - 1, 0) * factor * stride],
                                        line[c * factor * stride],
                                        line[(c + 1) * factor * stride],
                                        line[std::min<ptrdiff_t>(c + 2, n - 1) * factor * stride]};
                        for (ptrdiff_t m = 1; m < factor; ++m) {
                            const double *wm = &w[4 * m];
                            line[(c * factor + m) * stride] =
                                T(wm[0] * v[0] + wm[1] * v[1] + wm[2] * v[2] + wm[3] * v[3]);
                        }
                    }
                }
            }
        }
    }

    template <typename G, typename T>
    ptrdiff_t level_set_triangles(const GridGeometry<G> &grid,
                                  const T *const         f,
                                  const T                level,
                                  std::vector<G>        &x,
                                  std::vector<G>        &y,
                                  std::vector<G>        &z) {
        const ptrdiff_t        nk = grid.n[2] - 1;
        std::vector<ptrdiff_t> count(nk + 1, 0);

        // The triangles per layer of cells, then each layer's into its place.
#pragma omp parallel for schedule(dynamic, 1)
        for (ptrdiff_t k = 0; k < nk; ++k) {
            Cell      cell;
            double    tri[2][3][3];
            ptrdiff_t c = 0;
            for (ptrdiff_t j = 0; j + 1 < grid.n[1]; ++j) {
                for (ptrdiff_t i = 0; i + 1 < grid.n[0]; ++i) {
                    if (!load_cell(grid, f, level, i, j, k, cell)) continue;
                    for (int t = 0; t < 6; ++t) c += tet_triangles<G>(cell, kTets[t], tri);
                }
            }
            count[k + 1] = c;
        }
        for (ptrdiff_t k = 0; k < nk; ++k) count[k + 1] += count[k];

        const ptrdiff_t first = ptrdiff_t(x.size()) / 3, n_triangles = count[nk];
        x.resize(x.size() + 3 * n_triangles);
        y.resize(y.size() + 3 * n_triangles);
        z.resize(z.size() + 3 * n_triangles);

#pragma omp parallel for schedule(dynamic, 1)
        for (ptrdiff_t k = 0; k < nk; ++k) {
            Cell      cell;
            double    tri[2][3][3];
            ptrdiff_t next = first + count[k];
            for (ptrdiff_t j = 0; j + 1 < grid.n[1]; ++j) {
                for (ptrdiff_t i = 0; i + 1 < grid.n[0]; ++i) {
                    if (!load_cell(grid, f, level, i, j, k, cell)) continue;
                    for (int t = 0; t < 6; ++t) {
                        const int m = tet_triangles<G>(cell, kTets[t], tri);
                        for (int s = 0; s < m; ++s, ++next) {
                            for (int a = 0; a < 3; ++a) {
                                x[3 * next + a] = G(tri[s][a][0]);
                                y[3 * next + a] = G(tri[s][a][1]);
                                z[3 * next + a] = G(tri[s][a][2]);
                            }
                        }
                    }
                }
            }
        }
        return n_triangles;
    }

    template <typename G, typename T>
    int rebuild_sdf(const GridGeometry<G> &grid,
                    const T *const         f,
                    const T                level,
                    const T                band,
                    const ptrdiff_t        stride,
                    T *const               out) {
        const GridGeometry<G> og = grid.strided(stride);
        std::vector<G>        x, y, z;
        const ptrdiff_t       n_triangles = level_set_triangles(grid, f, level, x, y, z);

        std::vector<cuBQL::Triangle> triangles(n_triangles);
        std::vector<cuBQL::box3f>    boxes(n_triangles);
        for (ptrdiff_t t = 0; t < n_triangles; ++t) {
            triangles[t] = cuBQL::Triangle(cuBQL::vec3f(x[3 * t], y[3 * t], z[3 * t]),
                                           cuBQL::vec3f(x[3 * t + 1], y[3 * t + 1], z[3 * t + 1]),
                                           cuBQL::vec3f(x[3 * t + 2], y[3 * t + 2], z[3 * t + 2]));
            boxes[t]     = triangles[t].bounds();
        }
        cuBQL::bvh3f bvh;
        if (n_triangles > 0) cuBQL::cpuBuilder(bvh, boxes.data(), boxes.size(), cuBQL::BuildConfig());

        // Each sample's distance, searched within band only; a sample whose |f - level| is at
        // least band is at least band away (f is 1-Lipschitz: a distance, or a min or max of
        // distances) and is not searched.
        const ptrdiff_t sy = grid.n[0], sz = grid.n[0] * grid.n[1];
#pragma omp parallel for collapse(2) schedule(dynamic, 1)
        for (ptrdiff_t k = 0; k < og.n[2]; ++k) {
            for (ptrdiff_t j = 0; j < og.n[1]; ++j) {
                for (ptrdiff_t i = 0; i < og.n[0]; ++i) {
                    const ptrdiff_t p = i + og.n[0] * (j + og.n[1] * k);
                    const T         v = f[stride * (i + j * sy + k * sz)] - level;
                    T               d = band;
                    if (n_triangles > 0 && std::fabs(v) < band) {
                        const cuBQL::vec3f q(float(double(og.o[0]) + double(i) * double(og.h[0])),
                                             float(double(og.o[1]) + double(j) * double(og.h[1])),
                                             float(double(og.o[2]) + double(k) * double(og.h[2])));
                        cuBQL::triangles::CPAT cpat;
                        cpat.runQuery(triangles.data(), bvh, q, float(band));
                        if (cpat.triangleIdx >= 0) d = std::min(band, T(std::sqrt(cpat.sqrDist)));
                    }
                    out[p] = v > 0 ? d : (v < 0 ? -d : T(0));
                }
            }
        }
        if (n_triangles > 0) cuBQL::cpu::freeBVH(bvh);
        return 0;
    }

    template <typename G, typename T>
    int round_sdf(const GridGeometry<G> &grid,
                  const T *const         f,
                  const T                radius,
                  const T                band,
                  const bool             open,
                  const bool             close,
                  const ptrdiff_t        stride,
                  T *const               out,
                  T *const               f0) {
        const ptrdiff_t n = grid.size();
        std::vector<T>  own0(f0 ? 0 : n), a(open || close ? n : 0), b(open && close ? n : 0);
        T *const        r0 = f0 ? f0 : own0.data();

        auto shifted = [n](T *const v, const T s) {
#pragma omp parallel for
            for (ptrdiff_t p = 0; p < n; ++p) v[p] += s;
        };

        int rc = rebuild_sdf(grid, f, T(0), band, 1, r0);
        if (rc) return rc;
        const T *current = r0;
        if (open) {
            rc = rebuild_sdf(grid, current, radius, band, 1, a.data());
            if (rc) return rc;
            shifted(a.data(), radius);
            current = a.data();
        }
        if (close) {
            T *const next = open ? b.data() : a.data();
            rc            = rebuild_sdf(grid, current, -radius, band, 1, next);
            if (rc) return rc;
            shifted(next, -radius);
            current = next;
        }
        return rebuild_sdf(grid, current, T(0), band, stride, out);
    }

    template struct GridGeometry<float>;
    template struct GridGeometry<double>;

    template void refine_field<float, float>(const GridGeometry<float> &, const float *const, const ptrdiff_t, float *const);
    template void refine_field<double, double>(const GridGeometry<double> &,
                                               const double *const,
                                               const ptrdiff_t,
                                               double *const);

    template ptrdiff_t level_set_triangles<float, float>(const GridGeometry<float> &,
                                                         const float *const,
                                                         const float,
                                                         std::vector<float> &,
                                                         std::vector<float> &,
                                                         std::vector<float> &);
    template ptrdiff_t level_set_triangles<double, double>(const GridGeometry<double> &,
                                                           const double *const,
                                                           const double,
                                                           std::vector<double> &,
                                                           std::vector<double> &,
                                                           std::vector<double> &);

    template int rebuild_sdf<float, float>(const GridGeometry<float> &,
                                           const float *const,
                                           const float,
                                           const float,
                                           const ptrdiff_t,
                                           float *const);
    template int rebuild_sdf<double, double>(const GridGeometry<double> &,
                                             const double *const,
                                             const double,
                                             const double,
                                             const ptrdiff_t,
                                             double *const);

    template int round_sdf<float, float>(const GridGeometry<float> &,
                                         const float *const,
                                         const float,
                                         const float,
                                         const bool,
                                         const bool,
                                         const ptrdiff_t,
                                         float *const,
                                         float *const);
    template int round_sdf<double, double>(const GridGeometry<double> &,
                                           const double *const,
                                           const double,
                                           const double,
                                           const bool,
                                           const bool,
                                           const ptrdiff_t,
                                           double *const,
                                           double *const);

}  // namespace ssdf
