/*
Authors: 
Sourena Mosleh
Franz Zotter

Email-address: 
sourena.mosleh@student.kug.ac.at
zotter@iem.at

Institut fuer elektronische Musik (IEM)
Date and place: 15.07.2024, Graz
*/

#include "../fftw/fftw3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void hann(float* w, int len){
    float a=M_PI/(2*(len));
    for(int n=0; n<len; n++)
    {
        w[n]=cos(a*n)*cos(a*n);
        w[n]=sqrt(w[n]);
        
    }

}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void addArray(float * array1, float * array2, int len, float * result){

    for (int i=0; i<len; i++)
    result[i]=array1[i]+array2[i];
}
float squaredArray(float* array, int len)
{
    float result=0;
    for(int i=0; i<len; i++)
    result+=array[i]*array[i];

    return result;
}

float squared2DArray(float **x, int len1, int len2){
    float result=0;
    for (int i=0; i<len1;i++){
        result+=squaredArray(x[i],len2);
    }
    return result;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void resetArray(float *in, int len){
    for (int i=0; i<len;i++){
        in[i]=0;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset2DArray(float **in, int len1, int len2){
    for (int i=0; i<len1;i++){
	resetArray(in[i],len2);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset3DArray(float ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
	reset2DArray(in[i],len2,len3);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void unit_impulse(float *x, int len, int n){
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
void unitImpulse2D(float **x, int len1, int len2, int k, int n){
    reset2DArray(x,len1,len2);
    if ((k<len1)&&(n<len2)) {
       x[k][n]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void unitImpulse3D(float ***x, int len1, int len2, int len3, int m, int n, int o){
    reset3DArray(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)&&(o<len3)) {
         x[m][n][o]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void Rectangle2D(float **x, int len1, int len2, int k){
    reset2DArray(x,len1,len2);
    if ((k<len1)) {
        for (int i=0; i<len2;i++)
        x[k][i]=1;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void OneImpulse3D(float ***x, int len1, int len2, int len3, int m, int n){
    reset3DArray(x,len1,len2,len3);
    if ((m<len1)&&(n<len2)) {
        for (int i=0; i<len3;i++)
        x[m][n][i]=1;
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyComplexArray(fftwf_complex *src, fftwf_complex *dst, int len){
    for (int i=0; i<len;i++){
        dst[i][0]=src[i][0];
        dst[i][1]=src[i][1];
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyArray(float *src, float *dst, int len){
        for (int i=0; i<len;i++){
        dst[i]=src[i];
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void copy2DArray(float **src, float **dst, int channels,  int src_offset, int dst_offset, int len)
{
    for (int ch=0; ch<channels;ch++)
    {
        for (int num_sample=0; num_sample<len;num_sample++)
        dst[ch][dst_offset+num_sample]=src[ch][src_offset+num_sample];
    }
}

void copyExcerpt(float *src, float *dst, int src_offset, int dst_offset, int len)
{
    for (int n=0; n<len; n++)
    {
        dst[dst_offset+n]=src[src_offset+n];
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void copyArrayWithGain(float *src, float *dst, int len, float gain){
    for (int i=0; i<len;i++){
        dst[i]=src[i]/gain;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void resetComplexArray(fftwf_complex *in, int len){
        for (int i=0; i<len;i++){
            in[i][0]=0;
            in[i][1]=0;
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset2DComplexArray(fftwf_complex **in, int len1, int len2){
    for (int i=0; i<len1;i++){
	  resetComplexArray(in[i],len2);
    }
}


/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset3DComplexArray(fftwf_complex ***in, int len1, int len2, int len3){
    for (int i=0; i<len1;i++){
	  reset2DComplexArray(in[i],len2,len3);
    }
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset4DComplexArray(fftwf_complex ****in, int len1, int len2, int len3, int len4){
    for (int i=0; i<len1;i++){
	  reset3DComplexArray(in[i],len2,len3,len4);
    }
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void reset5DComplexArray(fftwf_complex *****in, int len1, int len2, int len3, int len4, int len5){
    for (int i=0; i<len1;i++){
	  reset4DComplexArray(in[i],len2,len3,len4,len5);
    }
}

/* COMPLEX MULTIPLY ACCUMULATE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void freq_mul_acc(fftwf_complex *array1, fftwf_complex *array2, fftwf_complex *result, int len){
    for (int i=0; i< len; i++)
    {
        result[i][0]+=array1[i][0]*array2[i][0]-array1[i][1]*array2[i][1];
        result[i][1]+=array1[i][1]*array2[i][0]+array1[i][0]*array2[i][1];
    }
}

/* PRINTING AND ERROR COMPUTATION TO DEBUG */
/*-----------------------------------------------------------------------------------------------------------------------------*/

void print_signal (float* in, int len){
     for(int i=0; i<len; i++)
        printf("%3.1f ", in[i]);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void print2DSignal (float** in, int len1, int len2){
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
void print_complex_signal (fftwf_complex* in, int len){
     for(int i=0; i<len; i++)
        printf("(%3.1f+i%3.1f) ", in[i][0], in[i][1]);
    printf("\n");
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
float squared_error(float *x,float *y, int len){
    float error=0;
    for (int i=0; i<len;i++){
        error+=(x[i]-y[i])*(x[i]-y[i]);
    }
    return error;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
float squared2DError(float **x,float **y, int len1, int len2){
    float error=0;
    for (int i=0; i<len1;i++){
        error+=squared_error(x[i],y[i],len2);
    }
    return error;
}

/*MEMORY MANAGEMENT FUNCTIONS*/
/*allocate memory*/
/*complex arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/

fftwf_complex* new1DComplexArray(int I)
{
    fftwf_complex* x= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * I);
    return x;
}

fftwf_complex** new2DComplexArray(int I,int J)
{
    fftwf_complex** x= (fftwf_complex**) fftwf_malloc(sizeof(fftwf_complex*) * I);
    for(int i=0; i<I; i++)
        x[i]= new1DComplexArray(J);
    return x;
}

fftwf_complex*** new3DComplexArray(int I,int J,int K)
{
    fftwf_complex*** x= (fftwf_complex***) fftwf_malloc(sizeof(fftwf_complex**) * I);
    for(int i=0; i<I; i++)
        x[i]= new2DComplexArray(J,K);
    return x;
}

fftwf_complex**** new4DComplexArray(int I,int J,int K,int L){
    fftwf_complex**** x= (fftwf_complex****) fftwf_malloc(sizeof(fftwf_complex***) * I);
    for(int i=0; i<I; i++)
        x[i]= new3DComplexArray(J,K,L);
    return x;
}

fftwf_complex***** new5DComplexArray(int I,int J,int K,int L, int M){
    fftwf_complex***** x= (fftwf_complex*****) fftwf_malloc(sizeof(fftwf_complex****) * I);
    for(int i=0; i<I; i++)
        x[i]= new4DComplexArray(J,K,L,M);
    return x;
}




/*float arrays*/
/*-----------------------------------------------------------------------------------------------------------------------------*/


float* new1DArray(int I)
{
    float* x= (float*) fftwf_malloc(sizeof(float) * I);
    return x;
}

float** new2DArray(int I,int J)
{
    float** x= (float**) fftwf_malloc(sizeof(float*) * I);
    for(int i=0; i<I; i++)
        x[i]= new1DArray(J);
    return x;
}

float*** new3DArray(int I,int J,int K)
{
    float*** x= (float***) fftwf_malloc(sizeof(float**) * I);
    for(int i=0; i<I; i++)
        x[i]= new2DArray(J,K);
    return x;
}

float**** new4DArray(int I,int J,int K,int L)
{
    float**** x= (float****) fftwf_malloc(sizeof(float***) * I);
    for(int i=0; i<I; i++)
        x[i]= new3DArray(J,K,L);
    return x;
}
/*free memory*/
/*-----------------------------------------------------------------------------------------------------------------------------*/
void free1DComplexArray(fftwf_complex* x)
{
    fftwf_free(x);
}

void free2DComplexArray(fftwf_complex** x, int I)
{
    for (int i=0; i<I; i++)
        free1DComplexArray(x[i]);
    fftwf_free(x);
}

void free3DComplexArray(fftwf_complex*** x, int I,int J)
{
    for (int i=0; i<I; i++)
        free2DComplexArray(x[i],J);
    fftwf_free(x);
}

void free4DComplexArray(fftwf_complex**** x, int I,int J,int K)
{
    for (int i=0; i<I; i++)
        free3DComplexArray(x[i],J,K);
    fftwf_free(x);
}
void free5DComplexArray(fftwf_complex***** x, int I,int J,int K,int L)
{
    for (int i=0; i<I; i++)
        free4DComplexArray(x[i],J,K,L);
    fftwf_free(x);
}



void free1DArray(float* x)
{
    fftwf_free(x);
}

void free2DArray(float** x, int I)
{
    for (int i=0; i<I; i++)
        free1DArray(x[i]);
    fftwf_free(x);
}

void free3DArray(float*** x, int I, int J)
{
    for (int i=0; i<I; i++)
        free2DArray(x[i],J);
    fftwf_free(x);
}
