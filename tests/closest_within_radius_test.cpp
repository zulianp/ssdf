#include "bvh/bvh.hpp"

#include <cassert>
#include <iostream>

using G = float;
using T = float;
using I = int;

void smoke_test() {
    const ptrdiff_t nspoints = 6;
    const G sx[nspoints] = {0.0f, 1.0f, 0.0f, 3.0f, 4.0f, 3.0f};
    const G sy[nspoints] = {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f};
    const G sz[nspoints] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    const ptrdiff_t nselements = 2;
    const I s0[nselements] = {0, 3};
    const I s1[nselements] = {1, 4};
    const I s2[nselements] = {2, 5};

    const ptrdiff_t npoints = 4;
    const G x[npoints] = {0.25f, 0.25f, 3.25f, 2.90f};
    const G y[npoints] = {0.25f, 0.25f, 0.25f, 0.20f};
    const G z[npoints] = {0.05f, 0.20f, 0.05f, 0.00f};

    const ptrdiff_t radius_stride = 1;
    const T radius_squared[npoints] = {0.01f, 0.01f, 0.01f, 10.00f};

    I outtri[npoints] = {-2, -2, -2, -2};
    const int result = ssdf::closest_within_radius_bvh<G, T, I>(
        npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, radius_stride, radius_squared, outtri);

    assert(result == 0);
    assert(outtri[0] == 0);
    assert(outtri[1] == -1);
    assert(outtri[2] == 1);
    assert(outtri[3] == 1);

    std::cout << "closest_within_radius_bvh tests passed\n";
}

void self_contact_test() {
    const ptrdiff_t nspoints = 4;
    const G sx[nspoints] = {0.0f, 1.0f, 0.0f, 0.2f};
    const G sy[nspoints] = {0.0f, 0.0f, 1.0f, 0.2f};
    const G sz[nspoints] = {0.0f, 0.0f, 0.0f, 0.0f};

    const ptrdiff_t nselements = 2;
    const I s0[nselements] = {0, 1};
    const I s1[nselements] = {1, 2};
    const I s2[nselements] = {2, 3};

    const ptrdiff_t npoints = 4;
    const G x[npoints] = {0.0f, 1.0f, 0.0f, 0.2f};
    const G y[npoints] = {0.0f, 0.0f, 1.0f, 0.2f};
    const G z[npoints] = {0.0f, 0.0f, 0.0f, 0.0f};

    const ptrdiff_t radius_stride = 1;
    const T radius_squared[npoints] = {1.0f, 1.0f, 1.0f, 1.0f};

    I outtri[npoints] = {-2, -2, -2, -2};
    const int result = ssdf::closest_within_radius_bvh<G, T, I>(npoints,
                                                                x,
                                                                y,
                                                                z,
                                                                nselements,
                                                                s0,
                                                                s1,
                                                                s2,
                                                                nspoints,
                                                                sx,
                                                                sy,
                                                                sz,
                                                                radius_stride,
                                                                radius_squared,
                                                                outtri,
                                                                true);

    assert(result == 0);
    assert(outtri[0] == 1);
    assert(outtri[1] == -1);
    assert(outtri[2] == -1);
    assert(outtri[3] == 0);
}

int main() {
    smoke_test();
    self_contact_test();
    return 0;
}
