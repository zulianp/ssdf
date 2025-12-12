#include "ssdf.hpp"
namespace ssdf {

    template int edf<float, float, int>(const ptrdiff_t npoints,
                                        const float *const SSDF_RESTRICT x,
                                        const float *const SSDF_RESTRICT y,
                                        const float *const SSDF_RESTRICT z,
                                        const ptrdiff_t nselements,
                                        const int *const SSDF_RESTRICT s0,
                                        const int *const SSDF_RESTRICT s1,
                                        const int *const SSDF_RESTRICT s2,
                                        const ptrdiff_t nspoints,
                                        const float *const SSDF_RESTRICT sx,
                                        const float *const SSDF_RESTRICT sy,
                                        const float *const SSDF_RESTRICT sz,
                                        float *const SSDF_RESTRICT out);

    template int edf<double, double, int>(const ptrdiff_t npoints,
                                          const double *const SSDF_RESTRICT x,
                                          const double *const SSDF_RESTRICT y,
                                          const double *const SSDF_RESTRICT z,
                                          const ptrdiff_t nselements,
                                          const int *const SSDF_RESTRICT s0,
                                          const int *const SSDF_RESTRICT s1,
                                          const int *const SSDF_RESTRICT s2,
                                          const ptrdiff_t nspoints,
                                          const double *const SSDF_RESTRICT sx,
                                          const double *const SSDF_RESTRICT sy,
                                          const double *const SSDF_RESTRICT sz,
                                          double *const SSDF_RESTRICT out);

}  // namespace ssdf

extern "C" {

int ssdf_edf_f(const ptrdiff_t npoints,
               const float *const SSDF_RESTRICT x,
               const float *const SSDF_RESTRICT y,
               const float *const SSDF_RESTRICT z,
               const ptrdiff_t nselements,
               const int *const SSDF_RESTRICT s0,
               const int *const SSDF_RESTRICT s1,
               const int *const SSDF_RESTRICT s2,
               const ptrdiff_t nspoints,
               const float *const SSDF_RESTRICT sx,
               const float *const SSDF_RESTRICT sy,
               const float *const SSDF_RESTRICT sz,
               float *const SSDF_RESTRICT out) {
    return ssdf::edf<float, float, int>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
}

int ssdf_edf_d(const ptrdiff_t npoints,
               const double *const SSDF_RESTRICT x,
               const double *const SSDF_RESTRICT y,
               const double *const SSDF_RESTRICT z,
               const ptrdiff_t nselements,
               const int *const SSDF_RESTRICT s0,
               const int *const SSDF_RESTRICT s1,
               const int *const SSDF_RESTRICT s2,
               const ptrdiff_t nspoints,
               const double *const SSDF_RESTRICT sx,
               const double *const SSDF_RESTRICT sy,
               const double *const SSDF_RESTRICT sz,
               double *const SSDF_RESTRICT out) {
    return ssdf::edf<double, double, int>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
}
}
