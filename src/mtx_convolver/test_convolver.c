/*
Test script/sketch board for
Uniformly Partitioned, Time-Variant,
Multichannel-Input-Mulichannel-Output Block Convolution
(and because signal processing folks like incomprehensible
 abbreviations: UPTVMIMOBC, yeah!)

useful for all kinds of real-time processes
virtual acoustic reality,
6dof spatial audio rendering/auralization,
convolution reverb,
Ambisonic decoding to headphones (binaural),
Ambisonic encoding of compact spherical microphone arrays
Ambisonic decoding to compact spherical loudspeaker arrays
Wave-Field Synthesis processing,
dynamic Binaural rendering,
Beamforming with loudspeaker or micrphone arrays,
adaptive signal processing,
frequency response equalization,
...

Franz Zotter
Hannes Pescoller
Sourena Mosleh

Email-address:
zotter@iem.at
h.pescoller@kug.ac.at
sourena.mosleh@student.kug.ac.at

Institute of Electronic Music and Acoustics (IEM)
University of Music and Performing Arts Graz
2024, 2025
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
  float eCounter = 0;
  float e;

  int M = 4;
  int L = 6;
  int P = 3;
  int Hann_len = L;
  int INPUT = 2;
  int OUTPUT = 5;

  conv_data *conv = initConvolution(L, P, Hann_len, INPUT, OUTPUT);

  // printf("hello!");
  in = new2DArray(conv->num_inputs, conv->L * M);
  inh = new3DArray(conv->num_outputs, conv->num_inputs, conv->L * conv->P);
  out = new2DArray(conv->num_outputs, conv->L * M);
  target = new2DArray(conv->num_outputs, conv->L * M);
  targettemp = new1DArray(conv->L * M);

  intemp = new2DArray(conv->num_inputs, conv->L);
  outtemp = new2DArray(conv->num_outputs, conv->L);

  /*reset in ,inh, out and target*/

  reset3DArray(inh, conv->num_outputs, conv->num_inputs, conv->L * M);
  reset2DArray(in, conv->num_inputs, conv->L * M);
  reset2DArray(intemp, conv->num_inputs, conv->L);
  reset2DArray(out, conv->num_outputs, conv->L * M);
  reset2DArray(target, conv->num_outputs, conv->L * M);
  reset2DArray(outtemp, conv->num_outputs, conv->L);

  for (int out_ch = 0; out_ch < conv->num_outputs; out_ch++) {
    for (int in_ch = 0; in_ch < conv->num_inputs; in_ch++) {
      for (int x_imp = 0; x_imp < conv->L * M; x_imp++) {
        unitImpulse2D(in, conv->num_inputs, conv->L * M, in_ch, x_imp);
        for (int h_imp = 0; h_imp < conv->L * conv->P; h_imp++) {
          unitImpulse3D(inh, conv->num_outputs, conv->num_inputs,
                        conv->L * conv->P, out_ch, in_ch, h_imp);
          unitImpulse2D(target, conv->num_outputs, conv->L * M, out_ch,
                        x_imp + h_imp);

          reset2DArray(conv->x_old, conv->num_inputs, conv->L);
          reset3DComplexArray(conv->xf, conv->num_inputs, conv->P, conv->L + 1);
          setImpulseResponse(conv, inh);
          printf("\nIR with Input channel %d, and Output channel %d=", in_ch,
                 out_ch);
          print_signal(inh[out_ch][in_ch], conv->P * conv->L);

          for (int m = 0; m < M; m++) {
            copy2DArray(in, intemp, conv->num_inputs, m * conv->L, 0, conv->L);
            conv_process(conv, intemp, outtemp);
            copy2DArray(outtemp, out, conv->num_outputs, 0, m * conv->L,
                        conv->L);
          }

          if ((e = fabsf(squared2DError(out, target, conv->num_outputs,
                                        conv->L))) > 0.3f) {
            printf("error i:%d,o:%d,xt:%d,ht:%d,e:%f\n", in_ch, out_ch, x_imp,
                   h_imp, e);
            eCounter++;
          } else {
            printf("no error i:%d,o:%d,xt:%d,ht:%d,e:%f\n", in_ch, out_ch,
                   x_imp, h_imp, e);
          }

          /*
              printf("\nOutput signal: ");
              print2DSignal(out,conv->num_outputs,conv->L*M);
              printf("\ntarget: ");
              print2DSignal(target,conv->num_outputs,conv->L*M);
              printf("\nInput signal: ");
              print2DSignal(in,conv->num_inputs,conv->L*M);
          */
        }
      }
    }
  }
  float errorRatio = eCounter / (conv->num_outputs * conv->num_inputs *
                                 conv->L * M * conv->L * conv->P);
  printf("\n%f errors from %d", eCounter,
         conv->num_outputs * conv->num_inputs * conv->L * M * conv->L *
             conv->P);
  printf("\nerror ratio= %f", errorRatio);
  /*for(int out_ch=0; out_ch<conv->num_outputs;out_ch++) {
      for(int in_ch=0; in_ch<conv->num_inputs;in_ch++)
          //for(int x_imp=0;x_imp<conv->L*M;x_imp++)
          {
              //int x_imp=0;
              //unitImpulse2D(in,conv->num_inputs,conv->L*M,in_ch,x_imp);
              Rectangle2D(in,conv->num_inputs,conv->L*M,in_ch);
              //for(int h_imp=0;h_imp<conv->L*conv->P;h_imp++) {
                  int h_imp=0;
                  //unitImpulse3D(inh,conv->num_outputs,conv->num_inputs,conv->L*conv->P,
                                  //out_ch, in_ch, h_imp);
                  //unitImpulse2D(target,conv->num_outputs,conv->L*M,
                                  //out_ch,x_imp+h_imp);
                  Rectangle2D(target,conv->num_outputs,conv->L*M,
                                  out_ch);


                  reset2DArray(conv->x_old,conv->num_inputs,conv->L);
                  reset3DComplexArray(conv->xf,conv->num_inputs,conv->P,conv->L+1);
                  int counter=0;
                  for (int m=0; m<M; m++)
                  {

                      if(counter%NUM_CF==0)
                      unitImpulse3D(inh,conv->num_outputs,conv->num_inputs,conv->L*conv->P,out_ch,
     in_ch, h_imp); else
                      reset3DArray(inh,conv->num_outputs,conv->num_inputs,conv->L*M);

                      setImpulseResponse(conv,inh);
                      //printf("convolver switch after set:%d
     \n",conv->convolver_switch); copy2DArray(in, intemp, conv->num_inputs,
     m*conv->L, 0, conv->L); conv_process(conv, intemp, outtemp);
                      //printf("convolver switch after process:%d
     \n",conv->convolver_switch); copy2DArray(outtemp, out, conv->num_outputs,
     0, m*conv->L, conv->L); counter++;
                  }


                  if
     ((e=fabsf(conv->num_outputs*squaredArray(conv->w,conv->L)-squared2DError(out,target,conv->num_outputs,conv->L)))>0.1f)
                  {
                      printf("error
     i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                      //print2DSignal(out,conv->num_outputs,conv->L*M);
                      //print2DSignal(target,conv->num_outputs,conv->L*M);
                  }
                  else
                  {
                      printf("no error
     i:%d,o:%d,xt:%d,ht:%d,e:%f\n",in_ch,out_ch,x_imp,h_imp,e);
                  }
                  */
  // printf("\nOutput signal: ");
  // print2DSignal(out,conv->num_outputs,conv->L*M);
  // printf("\nInput signal: ");
  // print2DSignal(in,conv->num_inputs,conv->L*M);
  // printf("\nThe Hann window: ");
  // print_signal(conv->w, conv->L);
  //}

  // }

  // }*/

  freeConvolution(conv);

  return 0;
}
