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
     * @param layers Input/Output array of layer indices. Inner layer is 0, outer layer is 1. Outer layer are marked by incremeneted values.
     * @param max_layers Maximum number of layers. If -1, the maximum number of layers is 255, 255 is used for unassigned layers.
     */
    template <typename I>
    int layers(const ptrdiff_t nelements,
               const int nxe,
               I** const SSDF_RESTRICT elements,
               const ptrdiff_t npoints,
               uint8_t* const SSDF_RESTRICT layers,
               const int max_layers = 255);
}  // namespace ssdf

#endif
