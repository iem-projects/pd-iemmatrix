#ifdef HAVE_FFTW
#include <fftw3.h>

#undef STUB
#define STUB(x) \
  void* iemmatrix_ ## x() { return x; }

STUB(fftw_malloc);
STUB(fftw_free);

STUB(fftw_plan_dft_c2r_1d);
STUB(fftw_plan_dft_r2c_1d);
STUB(fftw_execute);
STUB(fftw_destroy_plan);
#endif
