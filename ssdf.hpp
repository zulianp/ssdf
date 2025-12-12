#ifndef SSDF_HPP
#define SSDF_HPP

#include "cell_list.hpp"

#include <stddef.h>

#define SSDF_RESTRICT __restrict

namespace ssdf {

	// out is of size npoints
	template<typename G, typename T, typename I>
	int sdf(
		const ptrdiff_t npoints,
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
		T *const SSDF_RESTRICT out)
	{
		// 1) Compute surf aabbs
		// 2) Find closets distance on the x dimension
		//  - Sort surface w.r.t x
		//  - compute cum max
		//  - for each point binary search cum max to find closest element from the left use computed distance to conservatively discard candidates
		//  - Then check elements on the right  use computed distance to conservatively discard candidates
		// 3) Repeat for y and z
	}

} // ssdf

#endif //SSDF_HPP
