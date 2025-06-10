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

#ifdef USE_FFTWF
static const t_fftwf_malloc my_malloc = fftwf_malloc;
static const t_fftwf_free my_free = fftwf_free;
#else
static t_fftwf_malloc my_malloc = malloc;
static t_fftwf_free my_free = free;
#endif

int IEMCONVOLVE(array_set_fftwf_functions) (const t_fftwf_functions*funs) {
#ifndef USE_FFTWF
  if(funs && funs->malloc)
    my_malloc = funs->malloc;
  if(funs && funs->free)
    my_free = funs->free;
#endif
  return 1;
}


/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(cos2win) (float* w, int len){
    float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=cosf(a*n)*cosf(a*n);
    }
}
void IEMCONVOLVE(sin2win) (float* w, int len){
    float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=sinf(a*n)*sinf(a*n);
    }
}
void IEMCONVOLVE(coswin) (float* w, int len){
    float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=cosf(a*n);
    }
}
void IEMCONVOLVE(sinwin) (float* w, int len){
    float a=M_PI/(2*(len));
    for(int n=0; n<len; n++) {
        w[n]=sinf(a*n);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(addArray) (float * array1, float * array2, int len, float * result){
    for (int i=0; i<len; i++)
    result[i]=array1[i]+array2[i];
}
float IEMCONVOLVE(squaredArray) (float* array, int len)
{
    float result=0;
    for(int i=0; i<len; i++)
    result+=array[i]*array[i];

    return result;
}

float IEMCONVOLVE(squared2DArray) (float **x, int len1, int len2){
    float result=0;
    for (int i=0; i<len1;i++){
      result+=IEMCONVOLVE(squaredArray) (x[i],len2);
    }
    return result;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetArray) (float *in, int len){
    for (int i=0; i<len;i++){
        in[i]=0;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DArray) (float **in, int len1, int len2){
    for (int i=0; i<len1;i++){
	IEMCONVOLVE(resetArray)(in[i],len2);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DArray) (float ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
	IEMCONVOLVE(reset2DArray)(in[i],len2,len3);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unit_impulse) (float *x, int len, int n){
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
void IEMCONVOLVE(unitImpulse2D) (float **x, int len1, int len2, int k, int n){
    IEMCONVOLVE(reset2DArray)(x,len1,len2);
    if ((k<len1)&&(n<len2)) {
       x[k][n]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(unitImpulse3D) (float ***x, int len1, int len2, int len3, int m, int n, int o){
    IEMCONVOLVE(reset3DArray)(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)&&(o<len3)) {
         x[m][n][o]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(Rectangle2D) (float **x, int len1, int len2, int k){
    IEMCONVOLVE(reset2DArray)(x,len1,len2);
    if ((k<len1)) {
        for (int i=0; i<len2;i++)
        x[k][i]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(OneImpulse3D) (float ***x, int len1, int len2, int len3, int m, int n){
    IEMCONVOLVE(reset3DArray)(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)) {
        for (int i=0; i<len3;i++)
        x[m][n][i]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyComplexArray) (fftwf_complex *src, fftwf_complex *dst, int len){
    for (int i=0; i<len;i++){
        dst[i][0]=src[i][0];
        dst[i][1]=src[i][1];
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArray) (float *src, float *dst, int len){
        for (int i=0; i<len;i++){
        dst[i]=src[i];
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copy2DArray) (float **src, float **dst, int channels,  int src_offset, int dst_offset, int len)
{
    for (int ch=0; ch<channels;ch++)
    {
        for (int num_sample=0; num_sample<len;num_sample++)
        dst[ch][dst_offset+num_sample]=src[ch][src_offset+num_sample];
    }
}

void IEMCONVOLVE(copyExcerpt) (float *src, float *dst, int src_offset, int dst_offset, int len)
{
    for (int n=0; n<len; n++)
    {
        dst[dst_offset+n]=src[src_offset+n];
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(copyArrayWithGain) (float *src, float *dst, int len, float gain){
    for (int i=0; i<len;i++){
        dst[i]=src[i]/gain;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(resetComplexArray) (fftwf_complex *in, int len){
        for (int i=0; i<len;i++){
            in[i][0]=0;
            in[i][1]=0;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset2DComplexArray) (fftwf_complex **in, int len1, int len2){
    for (int i=0; i<len1;i++){
	  IEMCONVOLVE(resetComplexArray)(in[i],len2);
    }
}


/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset3DComplexArray) (fftwf_complex ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
	  IEMCONVOLVE(reset2DComplexArray)(in[i],len2,len3);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset4DComplexArray) (fftwf_complex ****in, int len1, int len2, int len3, int len4){
    for (int i=0; i<len1;i++){
	  IEMCONVOLVE(reset3DComplexArray)(in[i],len2,len3,len4);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(reset5DComplexArray) (fftwf_complex *****in, int len1, int len2, int len3, int len4, int len5){
    for (int i=0; i<len1;i++){
	  IEMCONVOLVE(reset4DComplexArray)(in[i],len2,len3,len4,len5);
    }
}

/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(freq_mul_acc) (fftwf_complex *array1, fftwf_complex *array2, fftwf_complex *result, int len){
    for (int i=0; i< len; i++)
    {
        result[i][0]+=array1[i][0]*array2[i][0]-array1[i][1]*array2[i][1];
        result[i][1]+=array1[i][1]*array2[i][0]+array1[i][0]*array2[i][1];
    }
}

/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/

void IEMCONVOLVE(print_signal ) (float* in, int len){
     for(int i=0; i<len; i++)
        printf("%3.1f ", in[i]);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(print2DSignal ) (float** in, int len1, int len2){
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
void IEMCONVOLVE(print_complex_signal ) (fftwf_complex* in, int len){
     for(int i=0; i<len; i++)
        printf("(%3.1f+i%3.1f) ", in[i][0], in[i][1]);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
float IEMCONVOLVE(squared_error) (float *x,float *y, int len){
    float error=0;
    for (int i=0; i<len;i++){
        error+=(x[i]-y[i])*(x[i]-y[i]);
    }
    return error;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
float IEMCONVOLVE(squared2DError) (float **x,float **y, int len1, int len2){
    float error=0;
    for (int i=0; i<len1;i++){
      error+=IEMCONVOLVE(squared_error) (x[i],y[i],len2);
    }
    return error;
}

/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/

fftwf_complex* IEMCONVOLVE(new1DComplexArray) (int I)
{
    if(!my_malloc)return 0;
    fftwf_complex* x= (fftwf_complex*) my_malloc(sizeof(fftwf_complex) * I);
    return x;
}

fftwf_complex** IEMCONVOLVE(new2DComplexArray) (int I,int J)
{
    if(!my_malloc)return 0;
    fftwf_complex** x= (fftwf_complex**) my_malloc(sizeof(fftwf_complex*) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new1DComplexArray) (J);
    return x;
}

fftwf_complex*** IEMCONVOLVE(new3DComplexArray) (int I,int J,int K)
{
    if(!my_malloc)return 0;
    fftwf_complex*** x= (fftwf_complex***) my_malloc(sizeof(fftwf_complex**) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new2DComplexArray) (J,K);
    return x;
}

fftwf_complex**** IEMCONVOLVE(new4DComplexArray) (int I,int J,int K,int L){
    if(!my_malloc)return 0;
    fftwf_complex**** x= (fftwf_complex****) my_malloc(sizeof(fftwf_complex***) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new3DComplexArray) (J,K,L);
    return x;
}

fftwf_complex***** IEMCONVOLVE(new5DComplexArray) (int I,int J,int K,int L, int M){
    if(!my_malloc)return 0;
    fftwf_complex***** x= (fftwf_complex*****) my_malloc(sizeof(fftwf_complex****) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new4DComplexArray) (J,K,L,M);
    return x;
}




/*float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/


float* IEMCONVOLVE(new1DArray) (int I)
{
    if(!my_malloc)return 0;
    float* x= (float*) my_malloc(sizeof(float) * I);
    return x;
}

float** IEMCONVOLVE(new2DArray) (int I,int J)
{
    if(!my_malloc)return 0;
    float** x= (float**) my_malloc(sizeof(float*) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new1DArray) (J);
    return x;
}

float*** IEMCONVOLVE(new3DArray) (int I,int J,int K)
{
    if(!my_malloc)return 0;
    float*** x= (float***) my_malloc(sizeof(float**) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new2DArray) (J,K);
    return x;
}

float**** IEMCONVOLVE(new4DArray) (int I,int J,int K,int L)
{
    if(!my_malloc)return 0;
    float**** x= (float****) my_malloc(sizeof(float***) * I);
    for(int i=0; i<I; i++)
        x[i]= IEMCONVOLVE(new3DArray) (J,K,L);
    return x;
}
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void IEMCONVOLVE(free1DComplexArray) (fftwf_complex* x)
{
    if(!my_free)return;
    my_free(x);
}

void IEMCONVOLVE(free2DComplexArray) (fftwf_complex** x, int I)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free1DComplexArray)(x[i]);
    if(!my_free)return;
    my_free(x);
}

void IEMCONVOLVE(free3DComplexArray) (fftwf_complex*** x, int I,int J)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free2DComplexArray)(x[i],J);
    if(!my_free)return;
    my_free(x);
}

void IEMCONVOLVE(free4DComplexArray) (fftwf_complex**** x, int I,int J,int K)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free3DComplexArray)(x[i],J,K);
    if(!my_free)return;
    my_free(x);
}
void IEMCONVOLVE(free5DComplexArray) (fftwf_complex***** x, int I,int J,int K,int L)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free4DComplexArray)(x[i],J,K,L);
    if(!my_free)return;
    my_free(x);
}



void IEMCONVOLVE(free1DArray) (float* x)
{
    if(!my_free)return;
    my_free(x);
}

void IEMCONVOLVE(free2DArray) (float** x, int I)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free1DArray)(x[i]);
    if(!my_free)return;
    my_free(x);
}

void IEMCONVOLVE(free3DArray) (float*** x, int I, int J)
{
    for (int i=0; i<I; i++)
        IEMCONVOLVE(free2DArray)(x[i],J);
    if(!my_free)return;
    my_free(x);
}
