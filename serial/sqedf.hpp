#ifndef SSDF_SQEDF_HPP
#define SSDF_SQEDF_HPP

#ifndef SSDF_RESTRICT
#ifndef _WIN32
#define SSDF_RESTRICT __restrict__
#else
#define SSDF_RESTRICT __restrict
#endif
#endif


namespace ssdf {

    /**
     * @brief Compute unsigned point-to-surface square distances.
     *
     * @tparam G Geometry type (float/double) for coordinates.
     * @tparam T Output distance type.
     * @tparam I Index type for triangle vertices.
     * @param npoints Number of query points.
     * @param x,y,z Arrays of size npoints with point coordinates (SoA).
     * @param nselements Number of surface triangles.
     * @param s0,s1,s2 Arrays of size nselements with vertex indices (triangle list).
     * @param nspoints Number of surface vertices.
     * @param sx,sy,sz Arrays of size nspoints with surface vertex coordinates (SoA).
     * @param out Output array of size npoints. On input, values are treated as
     *            current best square distances (initialize to large values if unused).
     * @return int 0 on success.
     *
     * Notes:
     * - Distances are unsigned (closest-point Euclidean distance).
     * - OpenMP is used when enabled at build time for parallelism.
     * - High-performance EDF using a 2D cell grid + sorted dimension
     */
    template <typename G, typename T, typename I>
    int sqedf(const ptrdiff_t npoints,
                     const G *const SSDF_RESTRICT x,
                     const G *const SSDF_RESTRICT y,
                     const G *const SSDF_RESTRICT z,
                     const ptrdiff_t nselements,
                     const I *const SSDF_RESTRICT s0,
                     const I *const SSDF_RESTRICT s1,
                     const I *const SSDF_RESTRICT s2,
                     const ptrdiff_t nspoints,
                     const G *const SSDF_RESTRICT sx,
                     const G *const SSDF_RESTRICT sy,
                     const G *const SSDF_RESTRICT sz,
                     T *const SSDF_RESTRICT out) {
        // SSDF_TIMER(edf_celllist);
        if (nselements == 0 || nspoints == 0) return 0;

        // Triangle AABBs
        std::vector<G> tminx(nselements), tminy(nselements), tminz(nselements);
        std::vector<G> tmaxx(nselements), tmaxy(nselements), tmaxz(nselements);

        G gminx = std::numeric_limits<G>::max();
        G gminy = std::numeric_limits<G>::max();
        G gminz = std::numeric_limits<G>::max();

        G gmaxx = std::numeric_limits<G>::lowest();
        G gmaxy = std::numeric_limits<G>::lowest();
        G gmaxz = std::numeric_limits<G>::lowest();

        G max_extent = G(0);

#pragma omp parallel
        {
            G lminx = std::numeric_limits<G>::max();
            G lminy = std::numeric_limits<G>::max();
            G lminz = std::numeric_limits<G>::max();

            G lmaxx = std::numeric_limits<G>::lowest();
            G lmaxy = std::numeric_limits<G>::lowest();
            G lmaxz = std::numeric_limits<G>::lowest();

            G lextent = G(0);
#pragma omp for nowait
            for (ptrdiff_t i = 0; i < nselements; ++i) {
                const I i0 = s0[i], i1 = s1[i], i2 = s2[i];
                const G x0 = sx[i0], x1 = sx[i1], x2 = sx[i2];
                const G y0 = sy[i0], y1 = sy[i1], y2 = sy[i2];
                const G z0 = sz[i0], z1 = sz[i1], z2 = sz[i2];

                const G xmin = std::min({x0, x1, x2});
                const G xmax = std::max({x0, x1, x2});
                const G ymin = std::min({y0, y1, y2});
                const G ymax = std::max({y0, y1, y2});
                const G zmin = std::min({z0, z1, z2});
                const G zmax = std::max({z0, z1, z2});

                tminx[i] = xmin;
                tmaxx[i] = xmax;
                tminy[i] = ymin;
                tmaxy[i] = ymax;
                tminz[i] = zmin;
                tmaxz[i] = zmax;

                lminx = std::min(lminx, xmin);
                lmaxx = std::max(lmaxx, xmax);
                lminy = std::min(lminy, ymin);
                lmaxy = std::max(lmaxy, ymax);
                lminz = std::min(lminz, zmin);
                lmaxz = std::max(lmaxz, zmax);
                lextent = std::max(lextent, std::max({xmax - xmin, ymax - ymin, zmax - zmin}));
            }
#pragma omp critical
            {
                gminx = std::min(gminx, lminx);
                gmaxx = std::max(gmaxx, lmaxx);
                gminy = std::min(gminy, lminy);
                gmaxy = std::max(gmaxy, lmaxy);
                gminz = std::min(gminz, lminz);
                gmaxz = std::max(gmaxz, lmaxz);
                max_extent = std::max(max_extent, lextent);
            }
        }

        // Choose grid axes: two widest spans
        const G span_x = gmaxx - gminx;
        const G span_y = gmaxy - gminy;
        const G span_z = gmaxz - gminz;
        int axis0 = 0, axis1 = 1, sort_axis = 2;
        {
            if (span_x >= span_y && span_x >= span_z) {
                axis0 = 0;
                if (span_y >= span_z) {
                    axis1 = 1;
                    sort_axis = 2;
                } else {
                    axis1 = 2;
                    sort_axis = 1;
                }
            } else if (span_y >= span_x && span_y >= span_z) {
                axis0 = 1;
                if (span_x >= span_z) {
                    axis1 = 0;
                    sort_axis = 2;
                } else {
                    axis1 = 2;
                    sort_axis = 0;
                }
            } else {
                axis0 = 2;
                if (span_x >= span_y) {
                    axis1 = 0;
                    sort_axis = 1;
                } else {
                    axis1 = 1;
                    sort_axis = 0;
                }
            }
        }

        const G cell_size = std::max(max_extent, G(1e-8));
        const G axis_min0 = (axis0 == 0) ? gminx : (axis0 == 1) ? gminy : gminz;
        const G axis_min1 = (axis1 == 0) ? gminx : (axis1 == 1) ? gminy : gminz;
        const G axis_span0 = (axis0 == 0) ? span_x : (axis0 == 1) ? span_y : span_z;
        const G axis_span1 = (axis1 == 0) ? span_x : (axis1 == 1) ? span_y : span_z;

        const I ncell0 = std::max(1, I(std::ceil(axis_span0 / cell_size)));
        const I ncell1 = std::max(1, I(std::ceil(axis_span1 / cell_size)));
        const ptrdiff_t ncells = ptrdiff_t(ncell0) * ptrdiff_t(ncell1);

        const G *tmin_axis0 = (axis0 == 0) ? tminx.data() : (axis0 == 1) ? tminy.data() : tminz.data();
        const G *tmin_axis1 = (axis1 == 0) ? tminx.data() : (axis1 == 1) ? tminy.data() : tminz.data();

        std::vector<I> cell_counts(ncells, 0);
        std::vector<G> cell_minx(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_miny(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_minz(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_maxx(ncells, std::numeric_limits<G>::lowest());
        std::vector<G> cell_maxy(ncells, std::numeric_limits<G>::lowest());
        std::vector<G> cell_maxz(ncells, std::numeric_limits<G>::lowest());

        auto cell_id = [&](G ax0, G ax1) -> ptrdiff_t {
            I ix = I((ax0 - axis_min0) / cell_size);
            I iy = I((ax1 - axis_min1) / cell_size);
            ix = std::max(0, std::min(ix, ncell0 - 1));
            iy = std::max(0, std::min(iy, ncell1 - 1));
            return ptrdiff_t(iy) * ncell0 + ix;
        };

        // Count and aggregate cell AABBs
        for (ptrdiff_t i = 0; i < nselements; ++i) {
            const ptrdiff_t cid = cell_id(tmin_axis0[i], tmin_axis1[i]);
            cell_counts[cid] += 1;
            cell_minx[cid] = std::min(cell_minx[cid], tminx[i]);
            cell_maxx[cid] = std::max(cell_maxx[cid], tmaxx[i]);
            cell_miny[cid] = std::min(cell_miny[cid], tminy[i]);
            cell_maxy[cid] = std::max(cell_maxy[cid], tmaxy[i]);
            cell_minz[cid] = std::min(cell_minz[cid], tminz[i]);
            cell_maxz[cid] = std::max(cell_maxz[cid], tmaxz[i]);
        }

        // Prefix sums
        std::vector<ptrdiff_t> cell_ptr(ncells + 1, 0);
        for (ptrdiff_t c = 0; c < ncells; ++c) cell_ptr[c + 1] = cell_ptr[c] + cell_counts[c];

        std::vector<I> cell_idx(nselements);
        std::vector<ptrdiff_t> fill_ptr = cell_ptr;
        for (ptrdiff_t i = 0; i < nselements; ++i) {
            const ptrdiff_t cid = cell_id(tmin_axis0[i], tmin_axis1[i]);
            cell_idx[fill_ptr[cid]++] = I(i);
        }

        const G *sorted_min = (sort_axis == 0) ? tminx.data() : (sort_axis == 1) ? tminy.data() : tminz.data();
        const G *sorted_max = (sort_axis == 0) ? tmaxx.data() : (sort_axis == 1) ? tmaxy.data() : tmaxz.data();

// Sort per cell on sort_axis
#pragma omp parallel for
        for (ptrdiff_t c = 0; c < ncells; ++c) {
            const ptrdiff_t begin = cell_ptr[c];
            const ptrdiff_t end = cell_ptr[c + 1];
            if (begin == end) continue;

            std::sort(cell_idx.begin() + begin, cell_idx.begin() + end, [&](I a, I b) {
                return sorted_min[a] < sorted_min[b];
            });
        }

        {  // Permute boxes to avoid aabb indirections during queries
            std::vector<G> temp(nselements);
            G *boxes[6] = {tminx.data(), tminy.data(), tminz.data(), tmaxx.data(), tmaxy.data(), tmaxz.data()};

            for (int d = 0; d < 6; d++) {
                memcpy(temp.data(), boxes[d], sizeof(G) * nselements);
                for (ptrdiff_t i = 0; i < nselements; i++) {
                    boxes[d][i] = temp[cell_idx[i]];
                }
            }
        }

        // build cumulative maxima per cell
        std::vector<G> cum_max(nselements, G(0));
        for (ptrdiff_t c = 0; c < ncells; ++c) {
            const ptrdiff_t begin = cell_ptr[c];
            const ptrdiff_t end = cell_ptr[c + 1];
            if (begin == end) continue;

            G acc = sorted_max[begin];
            for (ptrdiff_t i = begin; i < end; ++i) {
                acc = std::max(acc, sorted_max[i]);
                cum_max[i] = acc;
            }
        }

        {
            // SSDF_TIMER(queries);

            // Query points
#pragma omp parallel for
            for (ptrdiff_t p = 0; p < npoints; ++p) {
                const G px = x[p], py = y[p], pz = z[p];
                T best_sq = out[p];

                // TODO: change the brute force by
                // 1) finding cell containing the point
                // 2) process content of the cell
                // 3) Check sourrounding cells outword ring by ring
                // Make sure that the loop ends when no improvements are possible

                // Iterate cells with conservative culling
                for (ptrdiff_t cid = 0; cid < ncells; ++cid) {
                    if (cell_counts[cid] == 0 || !aabb_can_improve<T>(px,
                                                                      py,
                                                                      pz,
                                                                      best_sq,
                                                                      cell_minx[cid],
                                                                      cell_maxx[cid],
                                                                      cell_miny[cid],
                                                                      cell_maxy[cid],
                                                                      cell_minz[cid],
                                                                      cell_maxz[cid])) {
                        continue;
                    }

                    const ptrdiff_t begin = cell_ptr[cid];
                    const ptrdiff_t end = cell_ptr[cid + 1];
                    const G pcoord = (sort_axis == 0) ? px : (sort_axis == 1) ? py : pz;
                    const ptrdiff_t left = std::lower_bound(sorted_min + begin, sorted_min + end, pcoord) - sorted_min;

                    // Scan left
                    for (ptrdiff_t i = (left > begin) ? left - 1 : begin; i >= begin; --i) {
                        const G margin = pcoord - cum_max[i];
                        if (margin >= G(0) && margin * margin >= best_sq) break;

                        if (aabb_can_improve<T>(
                                px, py, pz, best_sq, tminx[i], tmaxx[i], tminy[i], tmaxy[i], tminz[i], tmaxz[i])) {
                            const I tid = cell_idx[i];
                            const I i0 = s0[tid], i1 = s1[tid], i2 = s2[tid];
                            const G dist_sq = point_triangle_dist_sq(
                                px, py, pz, sx[i0], sy[i0], sz[i0], sx[i1], sy[i1], sz[i1], sx[i2], sy[i2], sz[i2]);
                            best_sq = std::min(best_sq, dist_sq);
                        }
                    }

                    // Scan right
                    for (ptrdiff_t i = left; i < end; ++i) {
                        const G margin = sorted_min[i] - pcoord;
                        if (margin >= G(0) && margin * margin >= best_sq) break;

                        if (!aabb_can_improve<T>(
                                px, py, pz, best_sq, tminx[i], tmaxx[i], tminy[i], tmaxy[i], tminz[i], tmaxz[i])) {
                            continue;
                        }

                        const I tid = cell_idx[i];
                        const I i0 = s0[tid], i1 = s1[tid], i2 = s2[tid];
                        const G dist_sq = point_triangle_dist_sq(
                            px, py, pz, sx[i0], sy[i0], sz[i0], sx[i1], sy[i1], sz[i1], sx[i2], sy[i2], sz[i2]);
                        if (dist_sq < best_sq) best_sq = dist_sq;
                    }
                }

                out[p] = best_sq;
            }
        }

        return 0;
    }

}  // namespace ssdf

#endif  // SSDF_SQEDF_HPP
