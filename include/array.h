#include "fftw3.h"

/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void hann(float* w, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void addArray(float * array1, float * array2, int len, float * result);
float squaredArray(float* array, int len);
float squared2DArray(float **x, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void resetArray(float *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset2DArray(float **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset3DArray(float ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void unit_impulse(float *x, int len, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void unitImpulse2D(float **x, int len1, int len2, int k, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void unitImpulse3D(float ***x, int len1, int len2, int len3, int m, int n, int o);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void Rectangle2D(float **x, int len1, int len2, int k);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void OneImpulse3D(float ***x, int len1, int len2, int len3, int m, int n);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyComplexArray(fftwf_complex *src, fftwf_complex *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyArray(float *src, float *dst, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void copy2DArray(float **src, float **dst, int channels,  int src_offset, int dst_offset, int len);
void copyExcerpt(float *src, float *dst, int src_offset, int dst_offset, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyArrayWithGain(float *src, float *dst, int len, float gain);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void resetComplexArray(fftwf_complex *in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset2DComplexArray(fftwf_complex **in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset3DComplexArray(fftwf_complex ***in, int len1, int len2, int len3);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset4DComplexArray(fftwf_complex ****in, int len1, int len2, int len3, int len4);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset5DComplexArray(fftwf_complex *****in, int len1, int len2, int len3, int len4, int len5);
/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void freq_mul_acc(fftwf_complex *array1, fftwf_complex *array2, fftwf_complex *result, int len);
/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void print_signal (float* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void print2DSignal (float** in, int len1, int len2);
/*-----------------------------------------------------------------------------------------------------------------------------*/
void print_complex_signal (fftwf_complex* in, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
float squared_error(float *x,float *y, int len);
/*-----------------------------------------------------------------------------------------------------------------------------*/
float squared2DError(float **x,float **y, int len1, int len2);
/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
fftwf_complex* new1DComplexArray(int I);
fftwf_complex** new2DComplexArray(int I,int J);
fftwf_complex*** new3DComplexArray(int I,int J,int K);
fftwf_complex**** new4DComplexArray(int I,int J,int K,int L);
fftwf_complex***** new5DComplexArray(int I,int J,int K,int L, int M);
/*float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
float* new1DArray(int I);
float** new2DArray(int I,int J);
float*** new3DArray(int I,int J,int K);
float**** new4DArray(int I,int J,int K,int L);
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void free1DComplexArray(fftwf_complex* x);
void free2DComplexArray(fftwf_complex** x, int I);
void free3DComplexArray(fftwf_complex*** x, int I,int J);
void free4DComplexArray(fftwf_complex**** x, int I,int J,int K);
void free5DComplexArray(fftwf_complex***** x, int I,int J,int K,int L);
void free1DArray(float* x);
void free2DArray(float** x, int I);
void free3DArray(float*** x, int I, int J);
