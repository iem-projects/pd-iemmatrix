#ifndef _iemmatrix_stub_fftwf_h_
#define _iemmatrix_stub_fftwf_h_
#ifndef HAVE_FFTWF

/* from fftw3.h */
typedef struct _fftwf_plan_private_tag *fftwf_plan;
typedef float fftwf_complex[2];

#ifndef FFTW_ESTIMATE
# define FFTW_ESTIMATE (1U << 6)
#endif


#endif /* FFTWF */
#endif /* _iemmatrix_stub_fftwf_h_ */
