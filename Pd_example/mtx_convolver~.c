/*
Uniformly Partitioned, Time-Variant,
Multichannel-Input-Mulichannel-Output Block Convolution
(and because signal processing folks like incomprehensible
 abbreviations: UPTVMIMOBC, yeah!)
mtx_convolver~
for Pure-Data (with cross-faded outputs when updated)

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

#include "array.h"
#include "convolver.h"
#include "m_pd.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "iemmatrix.h"
#include <stdio.h>
#ifdef _WIN32
/* or should we use the security enhanced _snprintf_s() ?? */
# define snprintf _snprintf
#endif

# include <g_canvas.h>
# define IEMMATRIX_HAVE_G_CANVAS 1


static t_class *mtx_convolver_tilde_class;
static t_class *mtx_convolver_tilde_mclass;

typedef struct _mtx_convolver_tilde {
  t_object x_obj;
  int ins; // given/instantiated/mc number of inputs
  int outs; // given/instantiated number of inputs
  t_sample **inout_buffers;
  float ***h; // 3D impulse response array h, h_num_ins x h_num_outs x h_len from array.h
  int h_num_ins;
  int h_num_outs;
  int h_len;
  t_canvas *x_canvas; // to open array3 file
  conv_data *conv; // convolver structure from convolver.h
  float **conv_input_buffer; // input buffers to call convolver
  float **conv_output_buffer; // output buffers to call convolver
  t_inlet **x_inlet; // object's input ports
  t_outlet **x_outlet; // object's output ports
  int blocksize;
  int multichannel_mode;
  int set_ir_at_dsp_start; // retarded updated of input IR on dsp start
} t_mtx_convolver_tilde;

int ceildiv(int a, int b) {
  if (b == 0)
    return 0;
  int c = a / b;
  if (c * b < a)
    return c + 1;
  return c;
}

t_int *mtx_convolver_tilde_perform(t_int *w) {
  t_mtx_convolver_tilde *x = (t_mtx_convolver_tilde *)w[1];
  int available_inputs = (x->h_num_ins < x->ins) ? x->h_num_ins : x->ins;
  int available_outputs = (x->h_num_outs < x->outs) ? x->h_num_outs : x->outs;
  //post("[mtx_convolver~]: available ins=%d, outs=%d", available_inputs,available_outputs);
  if (x->conv) { // Operation: copy input signals, convolve, copy output signals
    for (int i = 0; i < available_inputs; i++) { // copy from inout pd buffers
      float *in = x->conv_input_buffer[i];
      t_sample *pd_in = x->inout_buffers[i];
      for (int n = 0; n < x->blocksize; n++) {
        in[n] = (float)pd_in[n];
      }
    }
    for (int i=available_inputs; i < x->h_num_ins; i++) { // zeropad of unavailable convolver inputs
      float *in = x->conv_input_buffer[i];
      for (int n = 0; n < x->blocksize; n++) {
        in[n] = (float)0;
      }
    }
    if (getNewIR(x->conv))
      post("[mtx_convolver~]: have to crossfade to a new IR");
    /*post("L=%d, P=%d, ins=%d, outs=%d",x->conv->L,x->conv->P,x->conv->num_inputs,x->conv->num_outputs);
    post("conv=%d, inbuf=%d, outpuf=%d",x->conv,x->conv_input_buffer, x->conv_output_buffer);*/
    conv_process(x->conv, x->conv_input_buffer, x->conv_output_buffer);
    for (int i = 0; i <available_outputs; i++) {
      float *out = x->conv_output_buffer[i];
      t_sample *pd_out = x->inout_buffers[i + x->ins];
      for (int n = 0; n < x->blocksize; n++) {
        pd_out[n] = (t_sample)out[n];
      }
    }
  } else { // DEFAULT: No convolver, deliver zero output signals.
    for (int i = 0; i < x->outs; i++) {
      t_sample *pd_out = x->inout_buffers[i + x->ins];
      for (int n = 0; n < x->blocksize; n++) {
        pd_out[n] = (t_sample)0;
      }
    }
  }
  return (w + 2);
}

void mtx_convolver_init(t_mtx_convolver_tilde *x) {
  if (x->conv)
    freeConvolution(x->conv);
  x->conv = 0;
  if ((x->blocksize > 0) && (x->h_num_ins > 0) && (x->h_num_outs > 0) && (x->h_len >0)) {
    int P = ceildiv(x->h_len, x->blocksize);
    x->conv = initConvolution(x->blocksize, P, x->blocksize, x->h_num_ins, x->h_num_outs);
  } 
}

int mtx_convolver_resize(t_mtx_convolver_tilde *x) {
  if (x->conv) {
    int P = ceildiv(x->h_len,x->blocksize);
    if ((x->conv->L==x->blocksize)&&(x->conv->P==P)&&
        (x->conv->num_inputs==x->h_num_ins)&&
        (x->conv->num_outputs==x->h_num_outs)) { 
        /*post("[mtx_convolver~]: keeping convolver with ins=%d, outs=%d, partitions=%d, blocksize=%d",
          x->conv->num_inputs,x->conv->num_outputs,x->conv->P,x->conv->L);*/
        return 0; // keep convolver, it's size is perfect, exit!
    } else { // resize required, resize buffers and free convolver
      if (x->conv_output_buffer) {
        free2DArray(x->conv_output_buffer, x->conv->num_outputs);
        x->conv_output_buffer=0;
      }
      if (x->conv_input_buffer) {
        free2DArray(x->conv_input_buffer, x->conv->num_inputs);
        x->conv_input_buffer=0;
      }
      freeConvolution(x->conv);
      x->conv = 0;
      if ((x->blocksize)&&(x->h_num_ins)&&(x->h_num_outs)) {
        x->conv_input_buffer = new2DArray(x->h_num_ins, x->blocksize);
        x->conv_output_buffer = new2DArray(x->h_num_outs, x->blocksize);
      }
    } 
 } // new convolver:
 if ((x->blocksize)&&(x->h_num_ins)&&(x->h_num_outs)&&(x->h_len)) {
  mtx_convolver_init(x);
  post("[mtx_convolver~]: re-instantiated convolver with ins=%d, outs=%d, partitions=%d, blocksize=%d",
      x->conv->num_inputs,x->conv->num_outputs,x->conv->P,x->conv->L);
  if (!x->conv_input_buffer)
    x->conv_input_buffer = new2DArray(x->h_num_ins, x->blocksize);
  if (!x->conv_output_buffer)
    x->conv_output_buffer = new2DArray(x->h_num_outs, x->blocksize);
 }
 return 1;
}

void mtx_convolver_tilde_dsp(t_mtx_convolver_tilde *x, t_signal **sp) {
  int ins = x->ins;
  int outs = x->outs;
  if (x->multichannel_mode) { 
    // override with num input signals in MC mode
    // override with num input signals in MC mode
    ins = sp[0]->s_nchans;
    post("[mtx_convolver~]: multichannel mode, in=%d",ins);
    // set number of outs in MC mode
    if (x->h_num_outs) {
        outs = x->h_num_outs;
        post("[mtx_convolver~]: multichannel mode, outs=%d", outs);
        signal_setmultiout(&sp[1], outs);
    }
    if ((ins+outs != x->ins+x->outs)){
      if (x->inout_buffers) {
        free(x->inout_buffers);
        x->inout_buffers=0;
      }
    }
    x->outs = outs;
    x->ins = ins;
  } else {
    /* post("[mtx_convolver~]: non-multichannel mode, in=%d",ins);
    post("[mtx_convolver~]: non-multichannel mode, outs=%d", outs); */
  }
  if (!x->inout_buffers) {
      x->inout_buffers = (t_float **)malloc(sizeof(float *) * (x->ins + x->outs));
  }
  x->blocksize = sp[0]->s_n;
  //post("[mtx_convolver~]: outs=%d, ins=%d, blocksize=%d",outs,ins,x->blocksize);
  int resized;
  resized = mtx_convolver_resize(x); // check if renewed convolver is required
  if (x->multichannel_mode) {
    for (int i = 0; i < x->ins; i++) {
      x->inout_buffers[i] = sp[0]->s_vec + i * x->blocksize;
    }
    for (int i = 0; i < x->outs; i++) {
      x->inout_buffers[x->ins + i] = sp[1]->s_vec + i * x->blocksize;
    }
  } else {
    for (int i = 0; i < x->ins + x->outs; i++) {
      x->inout_buffers[i] = sp[i]->s_vec;
    }
  }
  if ((x->set_ir_at_dsp_start) || (resized)) {
      if (x->conv) {
        /*post("[mtx_convolver~]: convolver has ins=%d, outs=%d, partitions=%d, blocksize=%d",
          x->conv->num_inputs,x->conv->num_outputs,x->conv->P,x->conv->L);*/
        setImpulseResponseZeroPad(x->conv, x->h, x->h_len);
        x->set_ir_at_dsp_start = 0;
      }
  }
  dsp_add(mtx_convolver_tilde_perform, 1, x);
}

void mtx_convolver_tilde_free(t_mtx_convolver_tilde *x) {
  if (x->inout_buffers)
    free(x->inout_buffers);
  if (x->h)
    free3DArray(x->h, x->h_num_outs, x->h_num_ins);
  if (x->conv)
    freeConvolution(x->conv);
  if (x->conv_input_buffer)
    free2DArray(x->conv_input_buffer, x->h_num_ins);
  if (x->conv_output_buffer)
    free2DArray(x->conv_output_buffer, x->h_num_outs);
}

const char *mtx_convolver_objname(void *obj) {
  t_object*x = obj;
  t_symbol*s = gensym("");
  if(x && x->te_binbuf) {
    char buf[MAXPDSTRING];
    t_symbol*objsym = atom_getsymbol(binbuf_getvec(x->te_binbuf));
    if(snprintf(buf, MAXPDSTRING, "[%s]: ", objsym->s_name) > 0) {
      buf[MAXPDSTRING-1] = 0;
      s = gensym(buf);
    }
  }
  return s->s_name;
}

int mtx_convolver_check(void*object, int argc, t_atom*argv, unsigned int tests) {
  t_object*x = (t_object*)object;
  const char*objname=mtx_convolver_objname(x);
  int inputs=(argc>2)?atom_getfloat(argv+0):0;
  int outputs=(argc>2)?atom_getfloat(argv+1):0;
  int length=(argc>2)?atom_getfloat(argv+2):0;
  if (!tests)
    tests =
      IEMMATRIX_CHECK_CRIPPLED
      | IEMMATRIX_CHECK_DIMENSIONS
      | IEMMATRIX_CHECK_SPARSE;
  if ((tests & IEMMATRIX_CHECK_CRIPPLED) && argc<3) {
    pd_error(x, "%scrippled array3", objname);
    return 1;
  }
  if ((tests & IEMMATRIX_CHECK_DIMENSIONS) && ((inputs<1)||(outputs<1)||(length<1))) {
    pd_error(x, "%sinvalid dimensions %dx%dx%d", objname, inputs, outputs, length);
    return 1;
  }
  if ((tests & IEMMATRIX_CHECK_SPARSE)&&(inputs* outputs*length>argc-3)) {
    pd_error(x, "%ssparse array3 not yet supported", objname);
    return 1;
  }
  return 0;
}

void mtx_convolver_tilde_array3(t_mtx_convolver_tilde *x, t_symbol *s, int argc,
                      t_atom *argv) {
  int h_num_ins;
  int h_num_outs;
  int h_len;
  int resized_outs=0;
  if (argc < 3) {
    post("array3 message must have at least 3 arguments: num_inputs, num_outputs, ir_len");
    return;
  }
  if (mtx_convolver_check(x, argc, argv, 0)) return;
  h_num_ins = (int)atom_getfloat(argv);
  h_num_outs = (int)atom_getfloat(argv + 1);
  h_len = (int)atom_getfloat(argv + 2);
  argv += 3;
  //post("inputs=%d, outputs=%d, ir_len=%d, x->blocksize=%d", h_num_ins, h_num_outs, h_len, x->blocksize);

  if ((h_num_ins != x->h_num_ins) || (h_num_outs != x->h_num_outs) || (h_len != x->h_len)) {
    /*post("input array3: re-sizing/setting x->h_num_ins=%d, x->h_num_outs=%d, x->h_len=%d", h_num_ins,
         h_num_outs, h_len);  */
    if (x->h) 
      free3DArray(x->h, x->h_num_outs, x->h_num_ins);
    x->h = new3DArray(h_num_outs, h_num_ins, h_len);
    x->h_num_ins = h_num_ins;
    if (x->h_num_outs!=h_num_outs) {
      resized_outs=1;
    }
    x->h_num_outs = h_num_outs;
    x->h_len = h_len;
  } 
  for (int in = 0; in < x->h_num_ins; in++) { // store input array
    for (int out = 0; out < x->h_num_outs; out++) { 
      for (int i = 0; i < x->h_len; i++) {
        x->h[out][in][i] = (float)atom_getfloat(argv++);
      }
    }
  }
  if (x->blocksize) { // if blocksize is already known
    mtx_convolver_resize(x); // check if new convolver needs to be instantiated
    if (x->conv)  {
      if (canvas_dspstate) { // if convolver exists+dsp is running update IRs
        /*post("[mtx_convolver~]: attempting ir update, dsp is running");
        post("[mtx_convolver~]: convolver has ins=%d, outs=%d, partitions=%d, blocksize=%d",
          x->conv->num_inputs,x->conv->num_outputs,x->conv->P,x->conv->L);*/
        setImpulseResponseZeroPad(x->conv, x->h, x->h_len);
        if (resized_outs) {
          //post("[mtx_convolver~]: convolver changed, calling dsp startup");
          canvas_update_dsp();
        }
      } else {
        x->set_ir_at_dsp_start = 1;
        //post("[mtx_convolver~]: ir update postponed to dsp start, dsp is off");
      }
    } else {
        x->set_ir_at_dsp_start = 1;
        //post("[mtx_convolver~]: ir update postponed to dsp start, convolver missing");
    }
  } else {
    x->set_ir_at_dsp_start = 1;
    //post("[mtx_convolver~]: ir update postponed to dsp start, blocksize missing");
  }
}

void mtx_convolver_tilde_read(t_mtx_convolver_tilde *x, t_symbol *filename)
{
  t_binbuf *bbuf = binbuf_new();
  t_atom *ap;
  int n;
  if (binbuf_read_via_path(bbuf, filename->s_name,
                           canvas_getdir(x->x_canvas)->s_name, 0)) {
    pd_error(x,"[mtx_convolver~]: failed to read '%s'", filename->s_name);
  }
  ap=binbuf_getvec(bbuf);
  n =binbuf_getnatom(bbuf)-1;
  if ((ap->a_type == A_SYMBOL) &&
      (!strcmp(ap->a_w.w_symbol->s_name,"array3")
       || !strcmp(ap->a_w.w_symbol->s_name,"#array3")) ) {
    mtx_convolver_tilde_array3(x, gensym("array3"), n, ap+1);
  }
  binbuf_free(bbuf);
}

void *mtx_convolver_tilde_new(t_symbol *s, int argc, t_atom *argv) {
  t_mtx_convolver_tilde *x;
  t_class *selected_class = mtx_convolver_tilde_class;

  if (argc > 0 && argv[0].a_type == A_SYMBOL &&
      strcmp(argv[0].a_w.w_symbol->s_name, "-m") == 0) {
    selected_class = mtx_convolver_tilde_mclass;
  }

  x = (t_mtx_convolver_tilde *)pd_new(selected_class);
  x->multichannel_mode = (selected_class == mtx_convolver_tilde_mclass);
  if (!x->multichannel_mode) {
    if (argc < 1) {
      x->ins = x->outs = 1;
    } else if (argc < 2) {
      x->ins = x->outs = (int)atom_getfloat(argv);
    } else {
      x->ins = (int)atom_getfloat(argv);
      x->outs = (int)atom_getfloat(argv + 1);
    }
    for (int i = 0; i < x->ins; i++) {
      inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    }
    for (int i = 0; i < x->outs; i++) {
      outlet_new(&x->x_obj, &s_signal);
    }
  }
  else {
    x->ins = x->outs = 1;
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
  }
  x->x_canvas = canvas_getcurrent();
  x->conv_input_buffer = 0;
  x->conv_output_buffer = 0;
  x->inout_buffers = 0;
  x->h = 0;
  x->conv = 0;
  x->blocksize = 0;
  x->set_ir_at_dsp_start = 0;

  if ((x->multichannel_mode) && (argc > 1 && argv[1].a_type == A_SYMBOL)) {
    t_symbol *matrix_file = argv[1].a_w.w_symbol;
    mtx_convolver_tilde_read(x, matrix_file);
  }
  if ((!x->multichannel_mode) && (argc > 2 && argv[2].a_type == A_SYMBOL)) {
    t_symbol *matrix_file = argv[2].a_w.w_symbol;
    mtx_convolver_tilde_read(x, matrix_file);
  }

  return (void *)x;
}

void mtx_convolver_tilde_setup(void) {
  mtx_convolver_tilde_class =
      class_new(gensym("mtx_convolver~"), (t_newmethod)mtx_convolver_tilde_new,
                (t_method)mtx_convolver_tilde_free,
                sizeof(t_mtx_convolver_tilde), CLASS_DEFAULT, A_GIMME, 0);

  mtx_convolver_tilde_mclass = class_new(
      gensym("mtx_convolver~ -m"), (t_newmethod)mtx_convolver_tilde_new,
      (t_method)mtx_convolver_tilde_free, sizeof(t_mtx_convolver_tilde),
      CLASS_MULTICHANNEL, A_GIMME, 0);

  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_array3,
                  gensym("array3"), A_GIMME, 0);
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_array3,
  gensym("#array3"), A_GIMME, 0);

  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_array3,
                  gensym("array3"), A_GIMME, 0);
  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_read,
                gensym("read"), A_SYMBOL, 0);
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_read,
                gensym("read"), A_SYMBOL, 0);
}
