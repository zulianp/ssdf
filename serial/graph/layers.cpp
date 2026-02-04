#include "layers.hpp"

#include <cassert>
#include <cstdio>

namespace ssdf {
    template <typename I, typename COUNT, typename ELINDEX>
    int create_n2e(const ptrdiff_t nelements,
                   const ptrdiff_t nnodes,
                   const int nnodesxelem,
                   I *SSDF_RESTRICT *const SSDF_RESTRICT elems,
                   COUNT **out_n2eptr,
                   ELINDEX **out_elindex) {
        COUNT *n2eptr = (COUNT *)malloc((nnodes + 1) * sizeof(COUNT));
        memset(n2eptr, 0, (nnodes + 1) * sizeof(COUNT));

        int *book_keeping = (int *)malloc((nnodes) * sizeof(int));
        memset(book_keeping, 0, (nnodes) * sizeof(int));

        for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            for (ptrdiff_t i = 0; i < nelements; ++i) {
                assert(elems[edof_i][i] < nnodes);
                assert(elems[edof_i][i] >= 0);

                ++n2eptr[elems[edof_i][i] + 1];
            }
        }

        for (ptrdiff_t i = 0; i < nnodes; ++i) {
            n2eptr[i + 1] += n2eptr[i];
        }

        ELINDEX *elindex = (ELINDEX *)malloc(n2eptr[nnodes] * sizeof(ELINDEX));

        for (int edof_i = 0; edof_i < nnodesxelem; ++edof_i) {
            for (ptrdiff_t i = 0; i < nelements; ++i) {
                ELINDEX node = elems[edof_i][i];
                assert(n2eptr[node] + book_keeping[node] < n2eptr[node + 1]);
                elindex[n2eptr[node] + book_keeping[node]++] = i;
            }
        }

        free(book_keeping);

        *out_n2eptr = n2eptr;
        *out_elindex = elindex;

        return 0;
    }

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
                         I *const SSDF_RESTRICT queue) {
        int cursor = *inout_cursor;
        int queue_size = *inout_queue_size;

        while (cursor < queue_size) {
            I current = queue[(cursor++) % queue_capacity];
            I ebegin = n2eptr[current];
            I eend = n2eptr[current + 1];

            for (I k = ebegin; k < eend; ++k) {
                I e = elindex[k];

                uint8_t min_layer = 255;
                for (int i = 0; i < nxe; ++i) {
                    I neighbor = elements[i][e];
                    min_layer = std::min(min_layer, layers[neighbor]);
                }

                if (min_layer + 1 == max_layers) {
                    continue;
                }

                for (int i = 0; i < nxe; ++i) {
                    I neighbor = elements[i][e];

                    if (layers[neighbor] <= min_layer + 1) {
                        continue;
                    }

                    queue[(queue_size++) % queue_capacity] = neighbor;
                    layers[neighbor] = min_layer + 1;
                }
            }
        }

        // Return the new cursor and queue size
        *inout_cursor = cursor;
        *inout_queue_size = queue_size;
        return 0;
    }

    template <typename I>
    int layers(const ptrdiff_t nelements,
               const int nxe,
               I *SSDF_RESTRICT *const SSDF_RESTRICT elements,
               const ptrdiff_t npoints,
               uint8_t *const SSDF_RESTRICT layers,
               const int max_layers) {
        ptrdiff_t *n2eptr = nullptr;
        I *elindex = nullptr;
        create_n2e<I, ptrdiff_t, I>(nelements, npoints, nxe, elements, &n2eptr, &elindex);

        I *queue = (I *)malloc(npoints * sizeof(I));
        int queue_size = 0;

        // Multi-source BFS initialization: every node with layer==0 is a seed.
        for (ptrdiff_t i = 0; i < npoints; ++i) {
            if (layers[i] == 0) {
                queue[queue_size++] = static_cast<I>(i);
            } else {
                // Make sure that the point is not assigned to a layer
                layers[i] = 255;
            }
        }

        if (queue_size == 0) {
            fprintf(stderr, "No seed found for layers (expect at least one node with layer==0)\n");
            free(queue);
            free(n2eptr);
            free(elindex);
            return 1;
        }

        int cursor = 0;

        layers_iterative(nxe, elements, layers, max_layers, n2eptr, elindex, &cursor, &queue_size, npoints, queue);

        free(queue);
        free(n2eptr);
        free(elindex);
        return 0;
    }

    template int create_n2e<int, ptrdiff_t, int>(const ptrdiff_t nelements,
                                                 const ptrdiff_t nnodes,
                                                 const int nnodesxelem,
                                                 int *SSDF_RESTRICT *const SSDF_RESTRICT elems,
                                                 ptrdiff_t **out_n2eptr,
                                                 int **out_elindex);

    template int layers<int>(const ptrdiff_t nelements,
                             const int nxe,
                             int *SSDF_RESTRICT *const SSDF_RESTRICT elems,
                             const ptrdiff_t npoints,
                             uint8_t *const SSDF_RESTRICT layers,
                             const int max_layers);

    template int layers_iterative<int>(const int nxe,
                                       int *SSDF_RESTRICT *const SSDF_RESTRICT elements,
                                       uint8_t *const SSDF_RESTRICT layers,
                                       const int max_layers,
                                       ptrdiff_t *const SSDF_RESTRICT n2eptr,
                                       int *const SSDF_RESTRICT elindex,
                                       int *inout_cursor,
                                       int *inout_queue_size,
                                       int queue_capacity,
                                       int *const SSDF_RESTRICT queue);
}  // namespace ssdf
