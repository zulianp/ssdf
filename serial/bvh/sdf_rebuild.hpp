#ifndef SSDF_SDF_REBUILD_HPP
#define SSDF_SDF_REBUILD_HPP

#include <stddef.h>

#include <vector>

namespace ssdf {

    /// A uniform grid: n[d] samples along axis d at o[d] + i h[d], x fastest.
    template <typename G>
    struct GridGeometry {
        ptrdiff_t n[3];
        G         o[3];
        G         h[3];

        ptrdiff_t size() const { return n[0] * n[1] * n[2]; }
        /// Every stride-th sample of this grid, from the first.
        GridGeometry strided(const ptrdiff_t stride) const;
        /// factor - 1 samples between neighbours: strided(factor) of the result is this grid.
        GridGeometry refined(const ptrdiff_t factor) const;
    };

    /// The triangles of the level set {f = level} of the field f on `grid`, by marching tetrahedra:
    /// each cell is split into the six tetrahedra about its diagonal from corner (0, 0, 0) to
    /// (1, 1, 1), so that neighbouring cells split their shared face alike (no cracks) and no case is
    /// ambiguous; the crossing on an edge by linear interpolation from its lower grid index (shared
    /// edges give the same vertex). A triangle soup: three vertices per triangle, appended to x, y,
    /// z. Returns the number of triangles.
    template <typename G, typename T>
    ptrdiff_t level_set_triangles(const GridGeometry<G> &grid,
                                  const T *const     f,
                                  const T            level,
                                  std::vector<G>    &x,
                                  std::vector<G>    &y,
                                  std::vector<G>    &z);

    /// out <- the signed distance to the level set {f = level} on every stride-th sample of `grid`
    /// (out has grid.strided(stride).size() entries): the magnitude the exact distance (to the
    /// float rounding of the closest-point queries) to the triangles of level_set_triangles, at most
    /// `band`; the sign that of f - level at the sample, zero on the level set. f must be
    /// 1-Lipschitz (a distance, or a min or max of distances): samples with |f - level| >= band
    /// are set to +-band without a search.
    template <typename G, typename T>
    int rebuild_sdf(const GridGeometry<G> &grid,
                    const T *const         f,
                    const T                level,
                    const T                band,
                    const ptrdiff_t        stride,
                    T *const               out);

    /// fine <- f on `grid` refined by `factor` (grid.refined(factor), the original samples kept): along
    /// each axis in turn, the cubic through the four nearest samples (Lagrange; the end samples
    /// repeated past the grid's ends).
    template <typename G, typename T>
    void refine_field(const GridGeometry<G> &grid, const T *const f, const ptrdiff_t factor, T *const fine);

    /// Rounds the edges of the region f > 0 by `radius` and returns the signed distance to the result
    /// (positive inside it), each step a rebuild_sdf, R(f, c), within `band`:
    ///   f0 = R(f, 0)               the distance to f's zero level set (a CSG min or max is none)
    ///   f1 = R(f0, r) + r          opening, erode then dilate by r: rounds the convex edges
    ///   f2 = R(f1, -r) - r         closing, dilate then erode by r: rounds the concave edges
    ///   out = R(f2, 0)             on every stride-th sample: exact on both sides
    /// The opening removes parts thinner than 2 r, the closing fills gaps narrower than 2 r. An edge
    /// of the eroded or dilated set (the intermediate level set of f1 or f2) is chamfered by its
    /// cells, so the rounded surface converges at first order where those edges were (second order
    /// elsewhere): process a grid finer than the output (stride > 1) to reduce that. band must exceed
    /// 3 r for the intermediate steps; f0 (may be null) receives R(f, 0) on the full grid.
    template <typename G, typename T>
    int round_sdf(const GridGeometry<G> &grid,
                  const T *const         f,
                  const T                radius,
                  const T                band,
                  const bool             open,
                  const bool             close,
                  const ptrdiff_t        stride,
                  T *const               out,
                  T *const               f0);

}  // namespace ssdf

#endif  // SSDF_SDF_REBUILD_HPP
