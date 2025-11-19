#ifndef _iemmatrix_stub_fftw_h_
#define _iemmatrix_stub_fftw_h_

/* from fftw3.h */
typedef struct _fftw_plan_private_tag *fftw_plan;
typedef float fftw_complex[2];

#ifndef FFTW_ESTIMATE
# define FFTW_ESTIMATE (1U << 6)
#endif

#endif /* _iemmatrix_stub_fftw_h_ */
