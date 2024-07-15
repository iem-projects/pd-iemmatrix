/*
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

#include "../include/array.h"
#include "../include/convolver.h"
#include <math.h>


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
    float eCounter=0;
    float e;

    int M=4;
    int L=6;
    int P=3;
    int Hann_len=L;
    int INPUT=2;
    int OUTPUT=5;

    
    conv_data* conv=initConvolution(L,P,Hann_len,INPUT,OUTPUT);
    

    //printf("hello!");
    in = new2DArray(conv->INPUT_Channel_Number,conv->L*M);
    inh = new3DArray(conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*conv->P);
    out = new2DArray(conv->OUTPUT_Channel_Number,conv->L*M);
    target = new2DArray(conv->OUTPUT_Channel_Number,conv->L*M);
    targettemp = new1DArray(conv->L*M);

    intemp = new2DArray(conv->INPUT_Channel_Number,conv->L);
    outtemp = new2DArray(conv->OUTPUT_Channel_Number,conv->L);

    
    /*reset in ,inh, out and target*/

    reset3DArray(inh,conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*M);  
    reset2DArray(in,conv->INPUT_Channel_Number,conv->L*M);
    reset2DArray(intemp,conv->INPUT_Channel_Number,conv->L);
    reset2DArray(out,conv->OUTPUT_Channel_Number,conv->L*M);
    reset2DArray(target,conv->OUTPUT_Channel_Number,conv->L*M);
    reset2DArray(outtemp,conv->OUTPUT_Channel_Number,conv->L);   


    for(int out_ch=0; out_ch<conv->OUTPUT_Channel_Number;out_ch++) 
    {
        for(int in_ch=0; in_ch<conv->INPUT_Channel_Number;in_ch++)
        { 
            for(int x_imp=0;x_imp<conv->L*M;x_imp++) 
            {
                unitImpulse2D(in,conv->INPUT_Channel_Number,conv->L*M,in_ch,x_imp);
                for(int h_imp=0;h_imp<conv->L*conv->P;h_imp++) 
                {
                    unitImpulse3D(inh,conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*conv->P, 
				    out_ch, in_ch, h_imp);
                    unitImpulse2D(target,conv->OUTPUT_Channel_Number,conv->L*M,
				    out_ch,x_imp+h_imp);

                    reset2DArray(conv->x_old,conv->INPUT_Channel_Number,conv->L);
                    reset3DComplexArray(conv->xf,conv->INPUT_Channel_Number,conv->P,conv->L+1);
                    setImpulseResponse(conv,inh);
                    printf("\nIR with Input channel %d, and Output channel %d=",in_ch, out_ch);
                    print_signal(inh[out_ch][in_ch],conv->P*conv->L);

                    for (int m=0; m<M; m++)
                    {
                        copy2DArray(in, intemp, conv->INPUT_Channel_Number, m*conv->L, 0, conv->L);
                        conv_process(conv, intemp, outtemp);
                        copy2DArray(outtemp, out, conv->OUTPUT_Channel_Number, 0, m*conv->L, conv->L);  
                    }

                    if ((e=fabsf(squared2DError(out,target,conv->OUTPUT_Channel_Number,conv->L)))>0.3f)
                    {
                        printf("error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                        eCounter++;
                    }
                    else
                    {
                        printf("no error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                    }


                    /*
                        printf("\nOutput signal: ");
                        print2DSignal(out,conv->OUTPUT_Channel_Number,conv->L*M);
                        printf("\ntarget: ");
                        print2DSignal(target,conv->OUTPUT_Channel_Number,conv->L*M);
                        printf("\nInput signal: ");
                        print2DSignal(in,conv->INPUT_Channel_Number,conv->L*M);
                    */ 

                }
            }
        }
    }
    float errorRatio=eCounter/(conv->OUTPUT_Channel_Number*conv->INPUT_Channel_Number*conv->L*M*conv->L*conv->P);
    printf("\n%f errors from %d",eCounter,conv->OUTPUT_Channel_Number*conv->INPUT_Channel_Number*conv->L*M*conv->L*conv->P);
    printf("\nerror ratio= %f", errorRatio);
    /*for(int out_ch=0; out_ch<conv->OUTPUT_Channel_Number;out_ch++) {
        for(int in_ch=0; in_ch<conv->INPUT_Channel_Number;in_ch++) 
            //for(int x_imp=0;x_imp<conv->L*M;x_imp++) 
            {
                //int x_imp=0;
                //unitImpulse2D(in,conv->INPUT_Channel_Number,conv->L*M,in_ch,x_imp);
                Rectangle2D(in,conv->INPUT_Channel_Number,conv->L*M,in_ch);
                //for(int h_imp=0;h_imp<conv->L*conv->P;h_imp++) {
                    int h_imp=0;
                    //unitImpulse3D(inh,conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*conv->P, 
				    //out_ch, in_ch, h_imp);
                    //unitImpulse2D(target,conv->OUTPUT_Channel_Number,conv->L*M,
				    //out_ch,x_imp+h_imp);
                    Rectangle2D(target,conv->OUTPUT_Channel_Number,conv->L*M,
				    out_ch);

                    
                    reset2DArray(conv->x_old,conv->INPUT_Channel_Number,conv->L);
                    reset3DComplexArray(conv->xf,conv->INPUT_Channel_Number,conv->P,conv->L+1);
                    int counter=0;
                    for (int m=0; m<M; m++)
                    {
                        
                        if(counter%NUM_CF==0)
                        unitImpulse3D(inh,conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*conv->P,out_ch, in_ch, h_imp);
                        else
                        reset3DArray(inh,conv->OUTPUT_Channel_Number,conv->INPUT_Channel_Number,conv->L*M);

                        setImpulseResponse(conv,inh);
                        //printf("convolver switch after set:%d \n",conv->convolver_switch);
                        copy2DArray(in, intemp, conv->INPUT_Channel_Number, m*conv->L, 0, conv->L);
                        conv_process(conv, intemp, outtemp);
                        //printf("convolver switch after process:%d \n",conv->convolver_switch);
                        copy2DArray(outtemp, out, conv->OUTPUT_Channel_Number, 0, m*conv->L, conv->L);  
                        counter++;
                    }

                    /*
                    if ((e=fabsf(conv->OUTPUT_Channel_Number*squaredArray(conv->w,conv->L)-squared2DError(out,target,conv->OUTPUT_Channel_Number,conv->L)))>0.1f)
                    {
                        printf("error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                        //print2DSignal(out,conv->OUTPUT_Channel_Number,conv->L*M);
                        //print2DSignal(target,conv->OUTPUT_Channel_Number,conv->L*M);
                    }
                    else
                    {
                        printf("no error i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                    }
                    */
                       // printf("\nOutput signal: ");
                        //print2DSignal(out,conv->OUTPUT_Channel_Number,conv->L*M);
                        //printf("\nInput signal: ");
                        //print2DSignal(in,conv->INPUT_Channel_Number,conv->L*M);
                        //printf("\nThe Hann window: ");
                        //print_signal(conv->w, conv->L);
                //}

            
            
       // }
        

   // }*/

    freeConvolution(conv); 
    
    return 0;
}

