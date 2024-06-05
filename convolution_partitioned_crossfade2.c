#include "fftw/fftw3.h"
#include <windows.h> 
#include <math.h> 
#define pi 3.1415
#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)
#define true 1
#define false 0



typedef struct Conv_data{
    int L;  // signal block length (fft length = 2L)
    int P;  // number of convolution partitions (length of h)
    int INPUT_Channel_Number; //number of input channels
    int OUTPUT_Channel_Number;
    int Hann_len;
    _Bool convolver_switch; //parameter for swiching bitween different buffers

    float *xtemp; // 2L fft-input time-domain signal x=[xprev, xcurr]  
    float *htemp; // 2L time-domain impulse response input h=[h, 0]
    float *y; // 2L time-domain result y=[ycircular, ylinear]
    float *y_cf; // 2L time-domain result y=[ycircular, ylinear] 
    float **x_old; //previous input data
    float *w; // hann-function
    int current_rb; // current_rb ring buffer position
    int current_cf;
    fftwf_complex ***xf; // L+1 positive-half DFT, partition input ring buffer
    fftwf_complex *****hf; // L+1 positive-half DFT, partition stack of h
    fftwf_complex *yftemp;  // L+1 positive-half DFT output buffer
    fftwf_complex *xftemp; // L+1 FFTW buffer for previous+current block
    fftwf_complex *hftemp; // L+1 FFTW buffer for response partition
    fftwf_plan fftplan_xtemp; // FFTW DFT plan for previous+current block
    fftwf_plan fftplan_htemp; // FFTW DFT plan for impulse response
    fftwf_plan ifftplan_y;    // FFTW iDFT plan for current output signal
    fftwf_plan ifftplan_y_cf;
    
} conv_data;

/* HELPER FUNCTIONS, GENERATION, COPYING, RESETTING */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void hann(float* w, int len){
    float a=pi/(2*(len));
    for(int n=0; n<len; n++)
    w[n]=cos(a*n)*cos(a*n);
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void addArray(float * array1, float * array2, int len, float * result){

    for (int i=0; i<len; i++)
    result[i]=array1[i]+array2[i];
}
float squaredArray(float* array, int len)
{
    float result;
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
    fftwf_complex***** x= (fftwf_complex*****) fftwf_malloc(sizeof(fftwf_complex) * I*J*K*L*M);
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

/* crossfade functions */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void crossfade(float *y,float* y_new, float* w, int len)
{
    for(int n=0; n<len; n++)
    y[n]=y[n]*w[n]+y_new[n]*(1-w[n]);
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
void setNewIR(conv_data *conv, _Bool status){
    conv->convolver_switch=status;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
_Bool getNewIR(conv_data *conv){
    _Bool result;
    result= conv->convolver_switch;
    return result;
}
/*-----------------------------------------------------------------------------------------------------------------------------*/
/* PARTITIONED CONVOLUTION CORE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
void conv_process(conv_data *conv, float **in, float **out){

    for (int in_ch=0; in_ch<conv->INPUT_Channel_Number; in_ch++)
    {
        copyArray(conv->x_old[in_ch],conv->xtemp,conv->L); // save old block
        copyArray(in[in_ch],&conv->xtemp[conv->L],conv->L); // combine with new block
        fftwf_execute(conv->fftplan_xtemp); // perform FFT
        copyComplexArray(conv->xftemp,conv->xf[in_ch][conv->current_rb],conv->L+1); // write into ring buffer
        copyArray(in[in_ch],conv->x_old[in_ch],conv->L); // store previous for next call
    }

    for (int out_ch=0; out_ch<conv->OUTPUT_Channel_Number; out_ch++)
    {
        resetComplexArray(conv->yftemp,conv->L+1); // zero output of partition accumulator
        for (int in_ch=0; in_ch<conv->INPUT_Channel_Number; in_ch++)
        {
            for (int p=0, px=conv->current_rb; p<conv->P; p++, px++)
            {
                freq_mul_acc(conv->xf[in_ch][px%conv->P],conv->hf[conv->current_cf][out_ch][in_ch][p],conv->yftemp,conv->L+1);
            }
        }
                
        fftwf_execute(conv->ifftplan_y); // perform accumulated partitions iFFT
        conv->convolver_switch=1;
        if(getNewIR(conv)==true)
        {
            conv->current_cf=(conv->current_cf+1)%NUM_CF;
            resetComplexArray(conv->yftemp,conv->L+1);
            for (int in_ch=0; in_ch<conv->INPUT_Channel_Number; in_ch++)
            {
                for (int p=0, px=conv->current_rb; p<conv->P; p++)
                {
                    freq_mul_acc(conv->xf[in_ch][px%conv->P],conv->hf[conv->current_cf][out_ch][in_ch][p],conv->yftemp,conv->L+1);
                    px++;
                }
            }
            fftwf_execute(conv->ifftplan_y_cf);
            crossfade(conv->y+conv->L,conv->y_cf+conv->L,conv->w,conv->L);
            setNewIR(conv,false);  //update status
        }
         copyArrayWithGain(&conv->y[conv->L],out[out_ch],conv->L,2*conv->L); //second half is result   
    }
    conv->current_rb=(conv->current_rb+conv->P-1)%conv->P; // decrease ring buffer position
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
conv_data initConvolution(int L, int P, int Hann_len, int in_ch, int out_ch){
    conv_data conv;
    conv.L=L;
    conv.P=P;
    conv.INPUT_Channel_Number=in_ch;
    conv.OUTPUT_Channel_Number=out_ch;
    conv.xtemp = new1DArray(conv.L*2);
    conv.htemp = new1DArray(conv.L*2);
    conv.y = new1DArray(conv.L*2);
    conv.y_cf = new1DArray(conv.L*2);
    conv.w = new1DArray(conv.L);
    conv.current_rb=0;
    conv.current_cf=1;
    conv.Hann_len=Hann_len;
    resetArray(conv.xtemp,conv.L*2);
    resetArray(conv.htemp,conv.L*2);
    resetArray(&conv.htemp[conv.L],conv.L); // zero pad
    resetArray(conv.w,conv.L);
    hann(conv.w,conv.Hann_len);
    conv.convolver_switch=0;


    conv.x_old=new2DArray(conv.INPUT_Channel_Number,conv.L);
    reset2DArray(conv.x_old,conv.INPUT_Channel_Number,conv.L);

    conv.xf=new3DComplexArray(conv.INPUT_Channel_Number, conv.P, conv.L+1);
    reset3DComplexArray(conv.xf,conv.INPUT_Channel_Number, conv.P, conv.L+1);
    conv.hf=new5DComplexArray(NUM_CF,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number, conv.P, conv.L+1);
    reset5DComplexArray(conv.hf,NUM_CF,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number, conv.P, conv.L+1);

    conv.xftemp = new1DComplexArray(conv.L+1); 
    conv.hftemp = new1DComplexArray(conv.L+1); 
    conv.yftemp = new1DComplexArray(conv.L+1); 

    conv.fftplan_xtemp = fftwf_plan_dft_r2c_1d(conv.L*2, conv.xtemp, conv.xftemp, FFTW_ESTIMATE);
    conv.fftplan_htemp = fftwf_plan_dft_r2c_1d(conv.L*2, conv.htemp, conv.hftemp, FFTW_ESTIMATE);
    conv.ifftplan_y = fftwf_plan_dft_c2r_1d(conv.L*2, conv.yftemp, conv.y, FFTW_ESTIMATE);
    conv.ifftplan_y_cf = fftwf_plan_dft_c2r_1d(conv.L*2, conv.yftemp, conv.y_cf, FFTW_ESTIMATE);
    return conv;
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void freeConvolution(conv_data conv){
    fftwf_destroy_plan(conv.fftplan_xtemp);
    fftwf_destroy_plan(conv.fftplan_htemp);
    fftwf_destroy_plan(conv.ifftplan_y);
    fftwf_destroy_plan(conv.ifftplan_y_cf);

    free3DComplexArray(conv.xf,conv.INPUT_Channel_Number,conv.P);
    free5DComplexArray(conv.hf,NUM_CF,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.P);
    free2DArray(conv.x_old,conv.INPUT_Channel_Number);

    fftwf_free(conv.yftemp);
    fftwf_free(conv.xtemp);
    fftwf_free(conv.y);
    fftwf_free(conv.w);
    fftwf_free(conv.y_cf);
    fftwf_free(conv.htemp);
    fftwf_free(conv.hftemp); 
    fftwf_free(conv.xftemp);
}

/*-----------------------------------------------------------------------------------------------------------------------------*/
void setImpulseResponse(conv_data *conv, float***inh){
    setNewIR(conv,true);
    conv->convolver_switch=1;

    for (int out_ch=0; out_ch<conv->OUTPUT_Channel_Number; out_ch++)
    {
        for (int in_ch=0; in_ch<conv->INPUT_Channel_Number; in_ch++)
        {
           for (int partition=0; partition<conv->P; partition++)
           {
                copyArray(&inh[out_ch][in_ch][partition*conv->L],conv->htemp,conv->L);
                fftwf_execute(conv->fftplan_htemp);
                copyComplexArray(conv->hftemp,conv->hf[(conv->current_cf+1)%NUM_CF][out_ch][in_ch][partition],conv->L+1);
            }
        }
    }
}
/* MAIN TESTING ROUTINE */
/*-----------------------------------------------------------------------------------------------------------------------------*/
int main() {

     
    float **in;
    float ***inh;
    float **out; 
    float **target;
    float *targettemp;
    float **intemp;
    float **outtemp;
    float e;

    int M=5;
    int L=6;
    int P=3;
    int Hann_len=L;
    int INPUT=2;
    int OUTPUT=5;

    
    conv_data conv=initConvolution(L,P,Hann_len,INPUT,OUTPUT);
    

    //printf("hello!");
    in = new2DArray(conv.INPUT_Channel_Number,conv.L*M);
    inh = new3DArray(conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.L*conv.P);
    out = new2DArray(conv.OUTPUT_Channel_Number,conv.L*M);
    target = new2DArray(conv.OUTPUT_Channel_Number,conv.L*M);
    targettemp = new1DArray(conv.L*M);

    intemp = new2DArray(conv.INPUT_Channel_Number,conv.L);
    outtemp = new2DArray(conv.OUTPUT_Channel_Number,conv.L);

    
    /*reset in ,inh, out and target*/

    reset3DArray(inh,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.L*M);  
    reset2DArray(in,conv.INPUT_Channel_Number,conv.L*M);
    reset2DArray(intemp,conv.INPUT_Channel_Number,conv.L);
    reset2DArray(out,conv.OUTPUT_Channel_Number,conv.L*M);
    reset2DArray(target,conv.OUTPUT_Channel_Number,conv.L*M);
    reset2DArray(outtemp,conv.OUTPUT_Channel_Number,conv.L);   

    
    for(int out_ch=0; out_ch<conv.OUTPUT_Channel_Number;out_ch++) {
        for(int in_ch=0; in_ch<conv.INPUT_Channel_Number;in_ch++) {
            //for(int x_imp=0;x_imp<conv.L*M;x_imp++) 
            {
                int x_imp=0;
                //unitImpulse2D(in,conv.INPUT_Channel_Number,conv.L*M,in_ch,x_imp);
                Rectangle2D(in,conv.INPUT_Channel_Number,conv.L*M,in_ch);
                //for(int h_imp=0;h_imp<conv.L*conv.P;h_imp++) {
                    int h_imp=0;
                    //unitImpulse3D(inh,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.L*conv.P, 
				    //out_ch, in_ch, h_imp);
                    //unitImpulse2D(target,conv.OUTPUT_Channel_Number,conv.L*M,
				    //out_ch,x_imp+h_imp);
                    Rectangle2D(target,conv.OUTPUT_Channel_Number,conv.L*M,
				    out_ch);

                    
                    reset2DArray(conv.x_old,conv.INPUT_Channel_Number,conv.L);
                    reset3DComplexArray(conv.xf,conv.INPUT_Channel_Number,conv.P,conv.L+1);
                    int counter=0;
                    for (int m=0; m<M; m++)
                    {
                        
                        if(counter%NUM_CF==0)
                        unitImpulse3D(inh,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.L*conv.P,out_ch, in_ch, h_imp);
                        else
                        reset3DArray(inh,conv.OUTPUT_Channel_Number,conv.INPUT_Channel_Number,conv.L*M);

                        setImpulseResponse(&conv,inh);
                        printf("convolver switch after set:%d \n",conv.convolver_switch);
                        copy2DArray(in, intemp, conv.INPUT_Channel_Number, m*conv.L, 0, conv.L);
                        conv_process(&conv, intemp, outtemp);
                        printf("convolver switch after process:%d \n",conv.convolver_switch);
                        copy2DArray(outtemp, out, conv.OUTPUT_Channel_Number, 0, m*conv.L, conv.L);  
                        counter++;
                    }
                    if ((e=fabsf(conv.OUTPUT_Channel_Number*squaredArray(conv.w,conv.L)-squared2DError(out,target,conv.OUTPUT_Channel_Number,conv.L)))>0.1f)
                    {
                        printf("error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                        //print2DSignal(out,conv.OUTPUT_Channel_Number,conv.L*M);
                        //print2DSignal(target,conv.OUTPUT_Channel_Number,conv.L*M);
                    }
                    else
                    {
                        printf("no error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                    }
                        print2DSignal(out,conv.OUTPUT_Channel_Number,conv.L*M);
                        print2DSignal(target,conv.OUTPUT_Channel_Number,conv.L*M);
                        print_signal(conv.w, conv.L);
                //}

            }
            
        }
        

    }

    freeConvolution(conv); 
    
    return 0;
}

