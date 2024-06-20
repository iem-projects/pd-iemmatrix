/*
 * HOWTO write an External for Pure data
 * (c) 2001-2006 IOhannes m zm�lnig zmoelnig[AT]iem.at
 *
 * this is the source-code for the fourth example in the HOWTO
 * it creates a simple dsp-object:
 * 2 input signals are mixed into 1 output signal
 * the mixing-factor can be set via the 3rd inlet
 *
 * for legal issues please see the file LICENSE.txt
 */


/**
 * include the interface to Pd
 */
#include "m_pd.h"
#include "convolver.h"
#include "array.h"
#include <stdlib.h>


/**
 * define a new "class"
 */
static t_class *mtx_convolver_tilde_class;


/**
 * this is the dataspace of our new object
 * the first element is the mandatory "t_object"
 * x_pan denotes the mixing-factor
 * "f" is a dummy and is used to be able to send floats AS signals.
 */
typedef struct _mtx_convolver_tilde {
  t_object x_obj;
  t_float x_pan;
  t_float f;
  int ins;
  int outs;
  int ins_old;
  int outs_old;
  float **hin; 
  float **input;
  float **output;
  int rows;
  int cols;
  conv_data *conv;
  t_inlet **x_in;
  t_outlet**x_out;
  int blocksize;
  int partition;
  t_signal **sp;
} t_mtx_convolver_tilde;


/**
 * this is the core of the object
 * this perform-routine is called for each signal block
 * the name of this function is arbitrary and is registered to Pd in the
 * mtx_convolver_tilde_dsp() function, each time the DSP is turned on
 *
 * the argument to this function is just a pointer within an array
 * we have to know for ourselves how many elements inthis array are
 * reserved for us (hint: we declare the number of used elements in the
 * mtx_convolver_tilde_dsp() at registration
 *
 * since all elements are of type "t_int" we have to cast them to whatever
 * we think is apropriate; "apropriate" is how we registered this function
 * in mtx_convolver_tilde_dsp()
 */
t_int *mtx_convolver_tilde_perform(t_int *w)
{
  t_mtx_convolver_tilde *x = (t_mtx_convolver_tilde*) (w+1); 
  t_signal **sp = (t_signal**) (w+2);

  // hier einfüllen: Signalblöcke dem laufenden 
  // Convolver (conv_process) übergeben
  // sp[0]->s_vec ist das t_sample (=float) Array vom Eingang 1
  // sp[1]->s_vec ist ... von Eingang 2...
  // muss herauskopiert werden.
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!out


  
  for(int i=0; i<x->ins;i++)
  {
    float* in=x->input[i];
    t_sample*pd_in = sp[i]->s_vec;
    for(int n=0; n<x->blocksize;n++)
    {
      in[n]=(float)pd_in[n];
    }
  }
  post("perf: x->conv=%d, x->input=%d. x->output=%d",x->conv,x->input,x->output);
/*
  conv_process(x->conv, x->input, x->output);

  for(int i=0; i<x->outs;i++)
  {
    float* out=x->output[i];
    t_sample*pd_out = sp[i+x->ins]->s_vec;
    for(int n=0; n<x->blocksize;n++)
    {
      pd_out[n]=(t_sample)out[n];
    }
  }
*/
 
  return (w+3);
}
void mtx_convolver_init(t_mtx_convolver_tilde *x,int ir_len)
{
  if (x->conv)
    freeConvolution(x->conv);
  if (x->blocksize) {
    int P = ceil(ir_len / x->blocksize);
    x->conv = initConvolution(x->blocksize, P, x->blocksize, x->ins, x->outs);
    post("init: x->conv=%d",x->conv);
  }
}


/**
 * register a special perform-routine at the dsp-engine
 * this function gets called whenever the DSP is turned ON
 * the name of this function is registered in mtx_convolver_tilde_setup()
 */
void mtx_convolver_tilde_dsp(t_mtx_convolver_tilde *x, t_signal **sp)
{
  /* add mtx_convolver_tilde_perform() to the DSP-tree;
   * the mtx_convolver_tilde_perform() will expect "5" arguments (packed into an
   * t_int-array), which are:
   * the objects data-space, 3 signal vectors (which happen to be
   * 2 input signals and 1 output signal) and the length of the
   * signal vectors (all vectors are of the same length)
   */

  // JMZ: here goes all the allocation (and pre-deallocation)

  if (x->output)
    free2DArray(x->output, x->outs);
  if (x->input)
    free2DArray(x->input, x->ins);

  x->input = new2DArray(x->ins, sp[0]->s_n);
  x->output = new2DArray(x->outs, sp[0]->s_n);

  post("dsp: x->cols=%d, x->blocksize=%d. sp[0]->s_n=%d", x->cols, x->blocksize, sp[0]->s_n);
  if ((x->blocksize != sp[0]->s_n) && (x->cols))
  {
    x->blocksize = sp[0]->s_n;
    int ir_len = ceil(x->cols / x->blocksize) * x->blocksize;
    mtx_convolver_init(x, ir_len);
    setImpulseResponse2DZeropad(x->conv, x->hin, x->cols);
    mtx_convolver_init(x, ir_len);

    post("dsp: x->conv=%d, x->input=%d. x->output=%d", x->conv, x->input, x->output);

    if (x->conv != 0)
      dsp_add(mtx_convolver_tilde_perform, 2, x, sp);
  }
  if (x->cols==0) {
    x->blocksize = sp[0]->s_n;
    x->sp=sp;
  }
}

/**
 * this is the "destructor" of the class;
 * it allows us to free dynamically allocated ressources
 */
void mtx_convolver_tilde_free(t_mtx_convolver_tilde *x)
{
  /* free any ressources associated with the given inlet */
  //inlet_free(x->x_in2);
  //inlet_free(x->x_in3);

  /* free any ressources associated with the given outlet */
 // outlet_free(x->x_out);
  if (x->input)
    free2DArray((float **)x->input, x->ins);
  if (x->output)
    free2DArray((float **)x->output, x->outs);
  if (x->hin)
     free2DArray(x->hin, x->rows);
  if (x->conv)
    freeConvolution(x->conv);
}

/**
 * this is the "constructor" of the class
 * the argument is the initial mixing-factor
 */
void *mtx_convolver_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
  t_mtx_convolver_tilde *x = (t_mtx_convolver_tilde *)pd_new(mtx_convolver_tilde_class);


  if (argc<1) {
    x->ins=x->outs=1;
  } else if (argc<2) {
    x->ins=x->outs=(int)atom_getfloat(argv);
  } else {
    x->ins=(int)atom_getfloat(argv);
    x->outs=(int)atom_getfloat(argv+1);
  }

 /*for (int i=0; i<argc; i++) {
    post("argv[%d]=%f",i,atom_getfloat(argv+i));
  }*/
  post("ins=%d, outs=%d",x->ins,x->outs);

  for(int i=0; i<x->ins; i++)
  {
     inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  }

  for(int i=0; i<x->outs; i++)
  {
     outlet_new(&x->x_obj, &s_signal);
  }
  x->input=0;
  x->output=0;
  x->hin=0;
  x->conv=0;
  x->blocksize=0;
  x->rows=0;
  x->cols=0;

  return (void *)x;
}


void mtx_input_matrix(t_mtx_convolver_tilde *x, t_symbol *s, int argc, t_atom *argv) {
  int rows, cols, ir_len=0;
  post("argc=%d",argc);
  //for(int i=0;i<argc;i++) {
   // post("%f ",atom_getfloat(argv+i));
 // }
  if (argc<3) {
     post("no valid matrix with just 2 entries");
     return;
  }
  rows=(int)atom_getfloat(argv);
  cols=(int)atom_getfloat(argv+1);
  if (rows*cols>argc-2) {
    post("no valid matrix: only %d elements instead of r x c=%d",argc-2,rows*cols);
     return;
  }
  if (rows!=x->ins*x->outs) {
    post("%d rows are not enough for %d inputs and %d outputs",rows,x->ins,x->outs);
    return;
  }
  if ((cols==0)||(rows==0)) {
    post("%d columns, %d rows specified", cols,rows);
    return;
  }
  argv+=2;
  //die Anzahl der Zeitsamples anpassen
  post("rows=%d, cols=%d, x->blocksize=%d",rows,cols,x->blocksize);
  if ((rows!=x->rows)||(cols!=x->cols)) {
    if (x->hin)
      free2DArray(x->hin,x->rows);
    x->hin=new2DArray(rows,cols);
    x->rows=rows;
    x->cols=cols;
  } 
  
  for (int io_idx = 0; io_idx < rows; io_idx++)
  {
      for (int time_index = 0; time_index < cols; time_index++)
      {
        x->hin[io_idx][time_index] = (float) atom_getfloat(argv++);
      }
  }
  if (x->blocksize)
  {
    ir_len = ceil(cols / x->blocksize) * x->blocksize;
    mtx_convolver_init(x, ir_len);
    setImpulseResponse2DZeropad(x->conv, x->hin, x->cols);
    dsp_add(mtx_convolver_tilde_perform, 2, x, sp);
  }
 // Registrieren eines Updates der Impulsantwort für den Convolver 
 // bzw: wenn Blockgröße bereits bekannt aus dsp-Routine
 // kann P usw. zum Initialisieren des Convolvers auch 
 // hier bereits gemacht werden, und setImpulseREsponse aufgerufen werden.
}


/**
 * define the function-space of the class
 * within a single-object external the name of this function is very special
 */
void mtx_convolver_tilde_setup(void) {
  mtx_convolver_tilde_class = class_new(gensym("mtx_convolver~"),
        (t_newmethod)mtx_convolver_tilde_new,
        (t_method)mtx_convolver_tilde_free,
        sizeof(t_mtx_convolver_tilde),
        CLASS_DEFAULT,
        A_GIMME, 0);

  /* whenever the audio-engine is turned on, the "mtx_convolver_tilde_dsp()"
   * function will get called
   */
  class_addmethod(mtx_convolver_tilde_class,
      (t_method)mtx_convolver_tilde_dsp, gensym("dsp"), A_CANT, 0);

  class_addmethod(mtx_convolver_tilde_class,
      (t_method)mtx_input_matrix, gensym("matrix"), 
      A_GIMME, 0);
  /* if no signal is connected to the first inlet, we can as well
   * connect a number box to it and use it as "signal"
   */
  CLASS_MAINSIGNALIN(mtx_convolver_tilde_class, t_mtx_convolver_tilde, f);
}
