#ifndef SSDF_BVH_HPP
#define SSDF_BVH_HPP

#include <stddef.h>

#include "ssdf_config.hpp"

#ifndef SSDF_RESTRICT
#ifndef _WIN32
#define SSDF_RESTRICT __restrict__
#else
#define SSDF_RESTRICT __restrict
#endif
#endif

namespace ssdf {

    template <typename G, typename T, typename I>
    int edf_bvh(const ptrdiff_t npoints,
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
                T *const SSDF_RESTRICT out);

}

#endif  // SSDF_BVH_HPP
