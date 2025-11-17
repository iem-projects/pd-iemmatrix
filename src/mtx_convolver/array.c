/*
Authors:
Sourena Mosleh
Franz Zotter

For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "LICENSE.txt," in this distribution.

Email-address:
sourena.mosleh@student.kug.ac.at
zotter@iem.at

Institute of Electronic Music and Acoustics (IEM)
University of Music and Performing Arts Graz
2024
*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif /* HAVE_CONFIG_H */

#include "array.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <iemmatrix_stub.h>
#include <iemmatrix_fft.h>


/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(cos2win) (t_float* w, int len){
    t_float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=cosf(a*n)*cosf(a*n);
    }
}
void IEMCONVOLVE(sin2win) (t_float* w, int len){
    t_float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=sinf(a*n)*sinf(a*n);
    }
}
void IEMCONVOLVE(coswin) (t_float* w, int len){
    t_float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=cosf(a*n);
    }
}
void IEMCONVOLVE(sinwin) (t_float* w, int len){
    t_float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=sinf(a*n);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(addArray) (const t_float * array1, const t_float * array2, int len, t_float * result){
    for (int i=0; i<len; i++)
    result[i]=array1[i]+array2[i];
}
t_float IEMCONVOLVE(squaredArray) (const t_float* array, int len)
{
    t_float result=0;
    for(int i=0; i<len; i++)
    result+=array[i]*array[i];

    return result;
}

t_float IEMCONVOLVE(squared2DArray) (const t_float **x, int len1, int len2){
    t_float result=0;
    for (int i=0; i<len1;i++){
      result+=IEMCONVOLVE(squaredArray) (x[i],len2);
    }
    return result;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetArray) (t_float *in, int len){
    for (int i=0; i<len;i++){
        in[i]=0;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DArray) (t_float **in, int len1, int len2){
    for (int i=0; i<len1;i++){
        IEMCONVOLVE(resetArray)(in[i],len2);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DArray) (t_float ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
        IEMCONVOLVE(reset2DArray)(in[i],len2,len3);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unit_impulse) (t_float *x, int len, int n){
    for (int i=0; i<len;i++){
        if (i==n){
            x[i]=1;
        }
        else{
            x[i]=0;
        }
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse2D) (t_float **x, int len1, int len2, int k, int n){
    IEMCONVOLVE(reset2DArray)(x,len1,len2);
    if ((k<len1)&&(n<len2)) {
       x[k][n]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse3D) (t_float ***x, int len1, int len2, int len3, int m, int n, int o){
    IEMCONVOLVE(reset3DArray)(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)&&(o<len3)) {
         x[m][n][o]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(Rectangle2D) (t_float **x, int len1, int len2, int k){
    IEMCONVOLVE(reset2DArray)(x,len1,len2);
    if ((k<len1)) {
        for (int i=0; i<len2;i++)
        x[k][i]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(OneImpulse3D) (t_float ***x, int len1, int len2, int len3, int m, int n){
    IEMCONVOLVE(reset3DArray)(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)) {
        for (int i=0; i<len3;i++)
        x[m][n][i]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyComplexArray) (const t_complex *src, t_complex *dst, int len){
    for (int i=0; i<len;i++){
        dst[i].re=src[i].re;
        dst[i].im=src[i].im;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArray) (const t_float *src, t_float *dst, int len){
        for (int i=0; i<len;i++){
        dst[i]=src[i];
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copy2DArray) (const t_float **src, t_float **dst, int channels,  int src_offset, int dst_offset, int len)
{
    for (int ch=0; ch<channels;ch++)
    {
        for (int num_sample=0; num_sample<len;num_sample++)
        dst[ch][dst_offset+num_sample]=src[ch][src_offset+num_sample];
    }
}

void IEMCONVOLVE(copyExcerpt) (const t_float *src, t_float *dst, int src_offset, int dst_offset, int len)
{
    for (int n=0; n<len; n++)
    {
        dst[dst_offset+n]=src[src_offset+n];
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArrayWithGain) (const t_float *src, t_float *dst, int len, t_float gain){
    for (int i=0; i<len;i++){
        dst[i]=src[i]/gain;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetComplexArray) (t_complex *in, int len){
        for (int i=0; i<len;i++){
            in[i].re=0;
            in[i].im=0;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DComplexArray) (t_complex **in, int len1, int len2){
    for (int i=0; i<len1;i++){
          IEMCONVOLVE(resetComplexArray)(in[i],len2);
    }
}


/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DComplexArray) (t_complex ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
          IEMCONVOLVE(reset2DComplexArray)(in[i],len2,len3);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset4DComplexArray) (t_complex ****in, int len1, int len2, int len3, int len4){
    for (int i=0; i<len1;i++){
          IEMCONVOLVE(reset3DComplexArray)(in[i],len2,len3,len4);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset5DComplexArray) (t_complex *****in, int len1, int len2, int len3, int len4, int len5){
    for (int i=0; i<len1;i++){
          IEMCONVOLVE(reset4DComplexArray)(in[i],len2,len3,len4,len5);
    }
}

/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freq_mul_acc) (const t_complex *array1, const t_complex *array2, t_complex *result, int len){
    for (int i=0; i< len; i++)
    {
        result[i].re+=array1[i].re*array2[i].re-array1[i].im*array2[i].im;
        result[i].im+=array1[i].im*array2[i].re+array1[i].re*array2[i].im;
    }
}

/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/

void IEMCONVOLVE(print_signal ) (const t_float* in, int len){
     for(int i=0; i<len; i++)
        printf("%3.1f ", in[i]);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print2DSignal ) (const t_float** in, int len1, int len2){
     for (int i=0; i<len1; i++) {
        if (i==0)
           printf("[");
        else
           printf(" ");
        for(int j=0; j<len2; j++) {
            printf("%3.1f ", in[i][j]);
        }
        if (i<len1-1)
           printf(";\n");
        else
           printf("]\n");
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print_complex_signal ) (const t_complex* in, int len){
     for(int i=0; i<len; i++)
        printf("(%3.1f+i%3.1f) ", in[i].re, in[i].im);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
t_float IEMCONVOLVE(squared_error) (const t_float *x, const t_float *y, int len){
    t_float error=0;
    for (int i=0; i<len;i++){
        error+=(x[i]-y[i])*(x[i]-y[i]);
    }
    return error;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
t_float IEMCONVOLVE(squared2DError) (const t_float **x, const t_float **y, int len1, int len2){
    t_float error=0;
    for (int i=0; i<len1;i++){
      error+=IEMCONVOLVE(squared_error) (x[i],y[i],len2);
    }
    return error;
}

/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/

t_complex* IEMCONVOLVE(new1DComplexArray) (int I)
{
    t_complex* x= (t_complex*) iemfft_malloc(sizeof(t_complex) * I);
    return x;
}

t_complex** IEMCONVOLVE(new2DComplexArray) (int I,int J)
{
    t_complex** x= (t_complex**) iemfft_malloc(sizeof(t_complex*) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new1DComplexArray) (J);
    return x;
}

t_complex*** IEMCONVOLVE(new3DComplexArray) (int I,int J,int K)
{
    t_complex*** x= (t_complex***) iemfft_malloc(sizeof(t_complex**) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new2DComplexArray) (J,K);
    return x;
}

t_complex**** IEMCONVOLVE(new4DComplexArray) (int I,int J,int K,int L){
    t_complex**** x= (t_complex****) iemfft_malloc(sizeof(t_complex***) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new3DComplexArray) (J,K,L);
    return x;
}

t_complex***** IEMCONVOLVE(new5DComplexArray) (int I,int J,int K,int L, int M){
    t_complex***** x= (t_complex*****) iemfft_malloc(sizeof(t_complex****) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new4DComplexArray) (J,K,L,M);
    return x;
}




/*t_float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/


t_float* IEMCONVOLVE(new1DArray) (int I)
{
    t_float* x= (t_float*) iemfft_malloc(sizeof(t_float) * I);
    return x;
}

t_float** IEMCONVOLVE(new2DArray) (int I,int J)
{
    t_float** x= (t_float**) iemfft_malloc(sizeof(t_float*) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new1DArray) (J);
    return x;
}

t_float*** IEMCONVOLVE(new3DArray) (int I,int J,int K)
{
    t_float*** x= (t_float***) iemfft_malloc(sizeof(t_float**) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new2DArray) (J,K);
    return x;
}

t_float**** IEMCONVOLVE(new4DArray) (int I,int J,int K,int L)
{
    t_float**** x= (t_float****) iemfft_malloc(sizeof(t_float***) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new3DArray) (J,K,L);
    return x;
}
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(free1DComplexArray) (t_complex* x)
{
    iemfft_free(x);
}

void IEMCONVOLVE(free2DComplexArray) (t_complex** x, int I)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free1DComplexArray)(x[i]);
    iemfft_free(x);
}

void IEMCONVOLVE(free3DComplexArray) (t_complex*** x, int I,int J)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free2DComplexArray)(x[i],J);
    iemfft_free(x);
}

void IEMCONVOLVE(free4DComplexArray) (t_complex**** x, int I,int J,int K)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free3DComplexArray)(x[i],J,K);
    iemfft_free(x);
}
void IEMCONVOLVE(free5DComplexArray) (t_complex***** x, int I,int J,int K,int L)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free4DComplexArray)(x[i],J,K,L);
    iemfft_free(x);
}



void IEMCONVOLVE(free1DArray) (t_float* x)
{
    iemfft_free(x);
}

void IEMCONVOLVE(free2DArray) (t_float** x, int I)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free1DArray)(x[i]);
    iemfft_free(x);
}

void IEMCONVOLVE(free3DArray) (t_float*** x, int I, int J)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free2DArray)(x[i],J);
    iemfft_free(x);
}
