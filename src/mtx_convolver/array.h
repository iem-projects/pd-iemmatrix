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

#include <m_pd.h>
#include <stddef.h>

#include <iemmatrix_fft.h>

#define IEMCONVOLVE(x) mtxconv_##x

/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(cos2win) (t_float* w, int len);
void IEMCONVOLVE(sin2win) (t_float* w, int len);
void IEMCONVOLVE(coswin) (t_float* w, int len);
void IEMCONVOLVE(sinwin) (t_float* w, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(addArray) (const t_float*array1, const t_float*array2, int len, t_float * result);
t_float IEMCONVOLVE(squaredArray) (const t_float* array, int len);
t_float IEMCONVOLVE(squared2DArray) (const t_float **x, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetArray) (t_float *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DArray) (t_float **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DArray) (t_float ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unit_impulse) (t_float *x, int len, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse2D) (t_float **x, int len1, int len2, int k, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse3D) (t_float ***x, int len1, int len2, int len3, int m, int n, int o);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(Rectangle2D) (t_float **x, int len1, int len2, int k);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(OneImpulse3D) (t_float ***x, int len1, int len2, int len3, int m, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyComplexArray) (const t_complex *src, t_complex *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArray) (const t_float *src, t_float *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copy2DArray) (const t_float **src, t_float **dst, int channels,  int src_offset, int dst_offset, int len);
void IEMCONVOLVE(copyExcerpt) (const t_float *src, t_float *dst, int src_offset, int dst_offset, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArrayWithGain) (const t_float *src, t_float *dst, int len, t_float gain);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetComplexArray) (t_complex *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DComplexArray) (t_complex **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DComplexArray) (t_complex ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset4DComplexArray) (t_complex ****in, int len1, int len2, int len3, int len4);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset5DComplexArray) (t_complex *****in, int len1, int len2, int len3, int len4, int len5);
/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freq_mul_acc) (const t_complex *array1, const t_complex *array2, t_complex *result, int len);
/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print_signal ) (const t_float* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print2DSignal ) (const t_float** in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print_complex_signal ) (const t_complex* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
t_float IEMCONVOLVE(squared_error) (const t_float *x, const t_float *y, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
t_float IEMCONVOLVE(squared2DError) (const t_float **x, const t_float **y, int len1, int len2);
/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
t_complex* IEMCONVOLVE(new1DComplexArray) (int I);
t_complex** IEMCONVOLVE(new2DComplexArray) (int I,int J);
t_complex*** IEMCONVOLVE(new3DComplexArray) (int I,int J,int K);
t_complex**** IEMCONVOLVE(new4DComplexArray) (int I,int J,int K,int L);
t_complex***** IEMCONVOLVE(new5DComplexArray) (int I,int J,int K,int L, int M);
/*t_float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
t_float* IEMCONVOLVE(new1DArray) (int I);
t_float** IEMCONVOLVE(new2DArray) (int I,int J);
t_float*** IEMCONVOLVE(new3DArray) (int I,int J,int K);
t_float**** IEMCONVOLVE(new4DArray) (int I,int J,int K,int L);
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(free1DComplexArray) (t_complex* x);
void IEMCONVOLVE(free2DComplexArray) (t_complex** x, int I);
void IEMCONVOLVE(free3DComplexArray) (t_complex*** x, int I,int J);
void IEMCONVOLVE(free4DComplexArray) (t_complex**** x, int I,int J,int K);
void IEMCONVOLVE(free5DComplexArray) (t_complex***** x, int I,int J,int K,int L);
void IEMCONVOLVE(free1DArray) (t_float* x);
void IEMCONVOLVE(free2DArray) (t_float** x, int I);
void IEMCONVOLVE(free3DArray) (t_float*** x, int I, int J);


#endif /* _mtx_convolver_array_h */
