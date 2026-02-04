#ifndef SSDF_GRAPH_LAYERS_HPP
#define SSDF_GRAPH_LAYERS_HPP

#include "ssdf.hpp"

namespace ssdf {
    /**
     * @brief Compute layers of a mesh.
     * @param nelements Number of elements.
     * @param nxe Number of nodes per element.
     * @param elements Arrays of element indices.
     * @param npoints Number of points.
     * @param layers Input/Output array of layer indices. Inner layer is 0, outer layer is 1. Outer layer are marked by
     * incremeneted values.
     * @param max_layers Maximum number of layers. If -1, the maximum number of layers is 255, 255 is used for
     * unassigned layers.
     */
    template <typename I>
    int layers(const ptrdiff_t nelements,
               const int nxe,
               I *SSDF_RESTRICT *const SSDF_RESTRICT elements,
               const ptrdiff_t npoints,
               uint8_t *const SSDF_RESTRICT layers,
               const int max_layers = 255);

    template <typename I, typename COUNT, typename ELINDEX>
    int create_n2e(const ptrdiff_t nelements,
                   const ptrdiff_t nnodes,
                   const int nnodesxelem,
                   I *SSDF_RESTRICT *const SSDF_RESTRICT elems,
                   COUNT **out_n2eptr,
                   ELINDEX **out_elindex);

    template <typename I>
    int layers_iterative(const int nxe,
                         I *SSDF_RESTRICT *const SSDF_RESTRICT elements,
                         uint8_t *const SSDF_RESTRICT layers,
                         const int max_layers,
                         //  Auxiliary arrays
                         ptrdiff_t *const SSDF_RESTRICT n2eptr,
                         I *const SSDF_RESTRICT elindex,
                         int *inout_cursor,
                         int *inout_queue_size,
                         int queue_capacity,
                         I *const SSDF_RESTRICT queue);
}  // namespace ssdf

#endif
