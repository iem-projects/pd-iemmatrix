/*
For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "LICENSE.txt," in this distribution.

Authors:
Sourena Mosleh
Franz Zotter

Email-address:
sourena.mosleh@student.kug.ac.at
zotter@iem.at

Institute of Electronic Music and Acoustics (IEM)
University of Music and Performing Arts Graz
2024
*/
#ifndef _mtx_convolver_array_h
#define _mtx_convolver_array_h

#include <stddef.h>

#if HAVE_FFTWF
# include <fftw3.h>
#else
# include "stub/fftwf.h"
#endif
typedef void*(*t_fftwf_malloc)(size_t);
typedef void(*t_fftwf_free)(void*);
typedef void(*t_fftwf_destroy_plan)(fftwf_plan);
typedef void(*t_fftwf_execute)(const fftwf_plan);
typedef fftwf_plan(*t_fftwf_plan_dft_r2c_1d)(int, float*, fftwf_complex*, unsigned);
typedef fftwf_plan(*t_fftwf_plan_dft_c2r_1d)(int, fftwf_complex*, float*, unsigned);

typedef struct _fftwf_functions {
  t_fftwf_malloc malloc;
  t_fftwf_free free;
  t_fftwf_destroy_plan destroy_plan;
  t_fftwf_execute execute;
  t_fftwf_plan_dft_r2c_1d plan_dft_r2c_1d;
  t_fftwf_plan_dft_c2r_1d plan_dft_c2r_1d;
} t_fftwf_functions;

#define IEMCONVOLVE(x) mtxconv_##x

/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
int IEMCONVOLVE(array_set_fftwf_functions) (const t_fftwf_functions*funs);

void IEMCONVOLVE(cos2win) (float* w, int len);
void IEMCONVOLVE(sin2win) (float* w, int len);
void IEMCONVOLVE(coswin) (float* w, int len);
void IEMCONVOLVE(sinwin) (float* w, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(addArray) (float * array1, float * array2, int len, float * result);
float IEMCONVOLVE(squaredArray) (float* array, int len);
float IEMCONVOLVE(squared2DArray) (float **x, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetArray) (float *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DArray) (float **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DArray) (float ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unit_impulse) (float *x, int len, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse2D) (float **x, int len1, int len2, int k, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse3D) (float ***x, int len1, int len2, int len3, int m, int n, int o);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(Rectangle2D) (float **x, int len1, int len2, int k);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(OneImpulse3D) (float ***x, int len1, int len2, int len3, int m, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyComplexArray) (fftwf_complex *src, fftwf_complex *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArray) (float *src, float *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copy2DArray) (float **src, float **dst, int channels,  int src_offset, int dst_offset, int len);
void IEMCONVOLVE(copyExcerpt) (float *src, float *dst, int src_offset, int dst_offset, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArrayWithGain) (float *src, float *dst, int len, float gain);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetComplexArray) (fftwf_complex *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DComplexArray) (fftwf_complex **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DComplexArray) (fftwf_complex ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset4DComplexArray) (fftwf_complex ****in, int len1, int len2, int len3, int len4);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset5DComplexArray) (fftwf_complex *****in, int len1, int len2, int len3, int len4, int len5);
/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freq_mul_acc) (fftwf_complex *array1, fftwf_complex *array2, fftwf_complex *result, int len);
/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print_signal ) (float* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print2DSignal ) (float** in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print_complex_signal ) (fftwf_complex* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
float IEMCONVOLVE(squared_error) (float *x,float *y, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
float IEMCONVOLVE(squared2DError) (float **x,float **y, int len1, int len2);
/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
fftwf_complex* IEMCONVOLVE(new1DComplexArray) (int I);
fftwf_complex** IEMCONVOLVE(new2DComplexArray) (int I,int J);
fftwf_complex*** IEMCONVOLVE(new3DComplexArray) (int I,int J,int K);
fftwf_complex**** IEMCONVOLVE(new4DComplexArray) (int I,int J,int K,int L);
fftwf_complex***** IEMCONVOLVE(new5DComplexArray) (int I,int J,int K,int L, int M);
/*float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
float* IEMCONVOLVE(new1DArray) (int I);
float** IEMCONVOLVE(new2DArray) (int I,int J);
float*** IEMCONVOLVE(new3DArray) (int I,int J,int K);
float**** IEMCONVOLVE(new4DArray) (int I,int J,int K,int L);
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(free1DComplexArray) (fftwf_complex* x);
void IEMCONVOLVE(free2DComplexArray) (fftwf_complex** x, int I);
void IEMCONVOLVE(free3DComplexArray) (fftwf_complex*** x, int I,int J);
void IEMCONVOLVE(free4DComplexArray) (fftwf_complex**** x, int I,int J,int K);
void IEMCONVOLVE(free5DComplexArray) (fftwf_complex***** x, int I,int J,int K,int L);
void IEMCONVOLVE(free1DArray) (float* x);
void IEMCONVOLVE(free2DArray) (float** x, int I);
void IEMCONVOLVE(free3DArray) (float*** x, int I, int J);


#endif /* _mtx_convolver_array_h */
