#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_bessel.h>

#undef STUB
#define STUB(x) \
  void *iemmatrix_##x() { return x; }

STUB(gsl_sf_bessel_Jn);
STUB(gsl_sf_bessel_Yn);

STUB(gsl_eigen_nonsymm);
STUB(gsl_eigen_nonsymm_alloc);
STUB(gsl_eigen_nonsymm_free);
STUB(gsl_eigen_nonsymmv);
STUB(gsl_eigen_nonsymmv_alloc);
STUB(gsl_eigen_nonsymmv_free);
STUB(gsl_linalg_QR_decomp);
STUB(gsl_linalg_QR_unpack);
STUB(gsl_linalg_SV_decomp);
STUB(gsl_matrix_alloc);
STUB(gsl_matrix_calloc);
STUB(gsl_matrix_set);
STUB(gsl_matrix_get);
STUB(gsl_matrix_complex_alloc);
STUB(gsl_matrix_complex_free);
STUB(gsl_matrix_free);
STUB(gsl_vector_alloc);
STUB(gsl_vector_complex_alloc);
STUB(gsl_vector_complex_free);
STUB(gsl_vector_free);
#endif
