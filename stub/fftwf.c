#ifdef HAVE_FFTWF
#include <fftw3.h>

#undef STUB
#define STUB(x) \
  void* iemmatrix_ ## x() { return x; }

STUB(fftwf_malloc);
STUB(fftwf_free);

STUB(fftwf_plan_dft_1d);
STUB(fftwf_plan_dft_c2r_1d);
STUB(fftwf_plan_dft_r2c_1d);
STUB(fftwf_execute);
STUB(fftwf_destroy_plan);
#endif
