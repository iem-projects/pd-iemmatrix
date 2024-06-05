#include "fftw/fftw3.h"
#include "array.h"
#define NUM_CF 2 // there are 2 crossfase buffers (re-occurring array dimension)



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

