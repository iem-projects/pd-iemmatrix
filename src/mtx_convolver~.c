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

For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "LICENSE.txt," in this distribution.

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

#include "iemmatrix.h"
#include "iemmatrix_stub.h"


#include "mtx_convolver/array.h"
#include "mtx_convolver/convolver.h"

#include <stdlib.h>

#if 1
static void _noppost(const void *object, int level, const char *fmt, ...) {
  (void)object; (void)level; (void)fmt;
}
# define _debug_logpost _noppost
#else
# define _debug_logpost logpost
#endif

#ifndef CLASS_MULTICHANNEL
# define CLASS_MULTICHANNEL 0
#endif

#if KERNEL_VERSION(PD_VERSION_MAJOR, PD_VERSION_MINOR, PD_VERSION_BUGFIX) < KERNEL_VERSION(0, 53, 0)
# define PD_NORMAL 2
#endif

typedef void (*setmultiout_f)(t_signal **sig, int nchans);

t_fftwf_functions my_functions;

static int warn_fftwf = 1;


static t_class *mtx_convolver_tilde_class;
static t_class *mtx_convolver_tilde_mclass;

typedef struct _mtx_convolver_tilde {
  t_object x_obj;
  t_symbol *x_objname;
  setmultiout_f x_setmultiout; /* when doing multichannel, this is Pd>=0.54's signal_setmultiout(); otherwise NULL */

  int ins; // given/instantiated/mc number of inputs
  int outs; // given/instantiated number of inputs
  t_sample **inout_buffers;
  float ***h; // 3D impulse response array h, h_num_ins x h_num_outs x h_len from array.h
  int h_num_ins;
  int h_num_outs;
  int h_len;
  _Bool coherent_xfade;
  t_canvas *x_canvas; // to open array3 file
  conv_data *conv; // convolver structure from convolver.h
  float **conv_input_buffer; // input buffers to call convolver
  float **conv_output_buffer; // output buffers to call convolver
  t_inlet **x_inlet; // object's input ports
  t_outlet **x_outlet; // object's output ports
  int blocksize;
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
  const char*objname=x->x_objname->s_name;
  int available_inputs = (x->h_num_ins < x->ins) ? x->h_num_ins : x->ins;
  int available_outputs = (x->h_num_outs < x->outs) ? x->h_num_outs : x->outs;
  _debug_logpost(x, PD_NORMAL, "[%s] available ins=%d, outs=%d", objname,available_inputs,available_outputs);
  if (x->conv) {
    _debug_logpost(x, PD_NORMAL, "[%s] processing with existing convolver",objname);
    // Operation: copy input signals, convolve, copy output signals
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
    if (IEMCONVOLVE(wasCrossFadeRegistered) (x->conv))
      logpost(x, PD_NORMAL, "[%s] have to crossfade to a new IR",objname);
    _debug_logpost(x, PD_NORMAL, "[%s] L=%d, P=%d, ins=%d, outs=%d",objname,x->conv->blocksize,x->conv->num_partitions,x->conv->num_inputs,x->conv->num_outputs);
    _debug_logpost(x, PD_NORMAL, "[%s] conv=%d, inbuf=%d, outpuf=%d",objname,x->conv,x->conv_input_buffer, x->conv_output_buffer);
    IEMCONVOLVE(convProcess) (x->conv, x->conv_input_buffer, x->conv_output_buffer);
    for (int i = 0; i <available_outputs; i++) {
      float *out = x->conv_output_buffer[i];
      t_sample *pd_out = x->inout_buffers[i + x->ins];
      for (int n = 0; n < x->blocksize; n++) {
        pd_out[n] = (t_sample)out[n];
      }
    }
  } else { // DEFAULT: No convolver, deliver zero output signals.
    _debug_logpost(x, PD_NORMAL, "[%s] NOP with convolver inexistent, zeroing %d outs with blocksize %d",objname,x->outs,x->blocksize);
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
    IEMCONVOLVE(freeConvolution) (x->conv);
  x->conv = 0;
  if ((x->blocksize > 0) && (x->h_num_ins > 0) && (x->h_num_outs > 0) && (x->h_len >0)) {
    int num_partitions = ceildiv(x->h_len, x->blocksize);
    x->conv = IEMCONVOLVE(initConvolution) (x->blocksize, num_partitions, x->blocksize, x->h_num_ins, x->h_num_outs, x->coherent_xfade);
  }
}

int mtx_convolver_resize(t_mtx_convolver_tilde *x) {
  const char*objname=x->x_objname->s_name;
  if (x->conv) {
    int num_partitions = ceildiv(x->h_len,x->blocksize);
    if ((x->conv->blocksize==x->blocksize)&&(x->conv->num_partitions==num_partitions)&&
        (x->conv->num_inputs==x->h_num_ins)&&
        (x->conv->num_outputs==x->h_num_outs)) {
      _debug_logpost(x, PD_NORMAL, "[%s] keeping convolver with ins=%d, outs=%d, partitions=%d, blocksize=%d",
			 objname,x->conv->num_inputs,x->conv->num_outputs,x->conv->num_partitions,x->conv->blocksize);
      return 0; // keep convolver, it's size is perfect, exit!
    } else { // resize required, resize buffers and free convolver
      if (x->conv_output_buffer) {
        IEMCONVOLVE(free2DArray) (x->conv_output_buffer, x->conv->num_outputs);
        x->conv_output_buffer=0;
      }
      if (x->conv_input_buffer) {
        IEMCONVOLVE(free2DArray) (x->conv_input_buffer, x->conv->num_inputs);
        x->conv_input_buffer=0;
      }
      IEMCONVOLVE(freeConvolution(x->conv));
      x->conv = 0;
      if ((x->blocksize)&&(x->h_num_ins)&&(x->h_num_outs)) {
        x->conv_input_buffer = IEMCONVOLVE(new2DArray) (x->h_num_ins, x->blocksize);
        x->conv_output_buffer = IEMCONVOLVE(new2DArray) (x->h_num_outs, x->blocksize);
      }
    }
 } // new convolver:
 if ((x->blocksize)&&(x->h_num_ins)&&(x->h_num_outs)&&(x->h_len)) {
  mtx_convolver_init(x);
  if (x->conv){
    logpost(x, PD_NORMAL, "[%s] re-instantiated convolver with ins=%d, outs=%d, partitions=%d, blocksize=%d",
      objname,x->conv->num_inputs,x->conv->num_outputs,x->conv->num_partitions,x->conv->blocksize);
  }
  if (!x->conv_input_buffer)
    x->conv_input_buffer = IEMCONVOLVE(new2DArray) (x->h_num_ins, x->blocksize);
  if (!x->conv_output_buffer)
    x->conv_output_buffer = IEMCONVOLVE(new2DArray) (x->h_num_outs, x->blocksize);
 }
 return 1;
}

void mtx_convolver_tilde_dsp(t_mtx_convolver_tilde *x, t_signal **sp) {
  int ins = x->ins;
  int outs = x->outs;
  const char*objname=x->x_objname->s_name;
#if CLASS_MULTICHANNEL
  if (x->x_setmultiout) {
    // override with num input signals in MC mode
    ins = sp[0]->s_nchans;
    logpost(x, PD_NORMAL, "[%s] multichannel mode, in=%d",objname,ins);
    // set number of outs in MC mode
    if (x->h_num_outs) {
        outs = x->h_num_outs;
        logpost(x, PD_NORMAL, "[%s] multichannel mode, outs=%d",objname,outs);
        x->x_setmultiout(&sp[1], outs);
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
    _debug_logpost(x, PD_NORMAL, "[%s] non-multichannel mode, in=%d, out=%d",objname,ins,outs);
  }
#endif
  if (!x->inout_buffers) {
      x->inout_buffers = (t_float **)malloc(sizeof(float *) * (x->ins + x->outs));
  }
  x->blocksize = sp[0]->s_n;
  _debug_logpost(x, PD_NORMAL, "[%s] outs=%d, ins=%d, blocksize=%d",objname,outs,ins,x->blocksize);
  int resized;
  resized = mtx_convolver_resize(x); // check if renewed convolver is required
  if (0) {
#if CLASS_MULTICHANNEL
  } else if (x->x_setmultiout) {
    for (int i = 0; i < x->ins; i++) {
      x->inout_buffers[i] = sp[0]->s_vec + i * x->blocksize;
    }
    for (int i = 0; i < x->outs; i++) {
      x->inout_buffers[x->ins + i] = sp[1]->s_vec + i * x->blocksize;
    }
#endif
  } else {
    for (int i = 0; i < x->ins + x->outs; i++) {
      x->inout_buffers[i] = sp[i]->s_vec;
    }
  }
  if ((x->set_ir_at_dsp_start) || (resized)) {
      if (x->conv) {
	_debug_logpost(x, PD_NORMAL, "[%s] convolver has ins=%d, outs=%d, partitions=%d, blocksize=%d",
			   objname,x->conv->num_inputs,x->conv->num_outputs,x->conv->num_partitions,x->conv->blocksize);
          IEMCONVOLVE(setImpulseResponseZeroPad) (x->conv, x->h, x->h_len, resized); // init with no xfade when conv was resized
        x->set_ir_at_dsp_start = 0;
      }
  }
  dsp_add(mtx_convolver_tilde_perform, 1, x);
}

void mtx_convolver_tilde_free(t_mtx_convolver_tilde *x) {
  if (x->inout_buffers)
    free(x->inout_buffers);
  if (x->h)
    IEMCONVOLVE(free3DArray) (x->h, x->h_num_outs, x->h_num_ins);
  if (x->conv)
    IEMCONVOLVE(freeConvolution) (x->conv);
  if (x->conv_input_buffer)
    IEMCONVOLVE(free2DArray) (x->conv_input_buffer, x->h_num_ins);
  if (x->conv_output_buffer)
    IEMCONVOLVE(free2DArray) (x->conv_output_buffer, x->h_num_outs);
}

int mtx_convolver_check(t_mtx_convolver_tilde*x, int argc, t_atom*argv, unsigned int tests) {
  const char*objname=x->x_objname->s_name;
  int inputs=(argc>2)?atom_getfloat(argv+0):0;
  int outputs=(argc>2)?atom_getfloat(argv+1):0;
  int length=(argc>2)?atom_getfloat(argv+2):0;
  if (!tests)
    tests =
      IEMMATRIX_CHECK_CRIPPLED
      | IEMMATRIX_CHECK_DIMENSIONS
      | IEMMATRIX_CHECK_SPARSE;
  if ((tests & IEMMATRIX_CHECK_CRIPPLED) && argc<3) {
    pd_error(x, "[%s] crippled array3", objname);
    return 1;
  }
  if ((tests & IEMMATRIX_CHECK_DIMENSIONS) && ((inputs<1)||(outputs<1)||(length<1))) {
    pd_error(x, "[%s] invalid dimensions %dx%dx%d", objname, inputs, outputs, length);
    return 1;
  }
  if ((tests & IEMMATRIX_CHECK_SPARSE)&&(inputs* outputs*length>argc-3)) {
    pd_error(x, "[%s] sparse array3 not yet supported", objname);
    return 1;
  }
  return 0;
}

void mtx_convolver_tilde_array3(t_mtx_convolver_tilde *x, t_symbol *s, int argc,
                      t_atom *argv) {
  const char*objname=x->x_objname->s_name;
  int h_num_ins;
  int h_num_outs;
  int h_len;
  int resized_outs=0;
  if (argc < 3) {
    logpost(x, PD_NORMAL, "[%s] %s message must have at least 3 arguments: num_inputs, num_outputs, ir_len",
            objname, s->s_name);
    return;
  }
  if (mtx_convolver_check(x, argc, argv, 0)) return;
  h_num_ins = (int)atom_getfloat(argv);
  h_num_outs = (int)atom_getfloat(argv + 1);
  h_len = (int)atom_getfloat(argv + 2);
  argv += 3;
  _debug_logpost(x, PD_NORMAL, "[%s] inputs=%d, outputs=%d, ir_len=%d, x->blocksize=%d", objname,h_num_ins, h_num_outs, h_len, x->blocksize);

  if ((h_num_ins != x->h_num_ins) || (h_num_outs != x->h_num_outs) || (h_len != x->h_len)) {
    _debug_logpost(x, PD_NORMAL, "[%s] input array3: re-sizing/setting x->h_num_ins=%d, x->h_num_outs=%d, x->h_len=%d", objname,h_num_ins,
		       h_num_outs, h_len);
    if (x->h)
      IEMCONVOLVE(free3DArray) (x->h, x->h_num_outs, x->h_num_ins);
    x->h = IEMCONVOLVE(new3DArray) (h_num_outs, h_num_ins, h_len);
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
    int resized = mtx_convolver_resize(x); // check if new convolver needs to be instantiated
    if (x->conv)  {
      if (canvas_dspstate) { // if convolver exists+dsp is running update IRs
	_debug_logpost(x, PD_NORMAL, "[%s] attempting IR update with dsp on, ins=%d, outs=%d, partitions=%d, blocksize=%d",
		       objname,x->conv->num_inputs,x->conv->num_outputs,x->conv->num_partitions,x->conv->blocksize);
        IEMCONVOLVE(setImpulseResponseZeroPad) (x->conv, x->h, x->h_len, resized); // init with no xfade if conv was resized, with xfade otherwise
        if (resized_outs) {
	  _debug_logpost(x, PD_NORMAL, "[%s] convolver changed, calling dsp startup",objname);
          canvas_update_dsp();
        }
      } else {
        x->set_ir_at_dsp_start = 1;
	_debug_logpost(x, PD_NORMAL, "[%s] ir update postponed to dsp start, dsp is off",objname);
      }
    } else {
      x->set_ir_at_dsp_start = 1;
      _debug_logpost(x, PD_NORMAL, "[%s] ir update postponed to dsp start, convolver missing",objname);
    }
  } else {
    x->set_ir_at_dsp_start = 1;
    _debug_logpost(x, PD_NORMAL, "[%s] ir update postponed to dsp start, blocksize missing",objname);
  }
}

void mtx_convolver_tilde_read(t_mtx_convolver_tilde *x, t_symbol *filename)
{
  const char *objname=x->x_objname->s_name;
  t_binbuf *bbuf = binbuf_new();
  t_atom *ap;
  int n;
  if (binbuf_read_via_path(bbuf, filename->s_name,
                           canvas_getdir(x->x_canvas)->s_name, 0)) {
    pd_error(x,"[%s] failed to read '%s'", objname, filename->s_name);
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
  setmultiout_f setmultiout = (CLASS_MULTICHANNEL)?iemmatrix_getpdfun("signal_setmultiout"):0;
  t_mtx_convolver_tilde *x = 0;
  t_class *selected_class = mtx_convolver_tilde_class;

  if (argc > 0 && argv[0].a_type == A_SYMBOL &&
      strcmp(argv[0].a_w.w_symbol->s_name, "-m") == 0) {
    /* want multichannel! */
#if CLASS_MULTICHANNEL
    if(setmultiout)
      selected_class = mtx_convolver_tilde_mclass;
    else {
      int major, minor, bugfix;
      sys_getversion(&major, &minor, &bugfix);
      pd_error(x, "[%s] multichannel requested, but iemmatrix is running in Pd-%d.%d-%d, which doesn't support it",
	       s->s_name, major, minor, bugfix);
      return 0;
    }

#else
    pd_error(x, "[%s] has been compiled without multichannel support", s->s_name);
    return 0;
#endif
  }

  x = (t_mtx_convolver_tilde *)pd_new(selected_class);
  x->x_objname = s;
  if(warn_fftwf) {
#if HAVE_FFTW
    pd_error(x, "[%s] couldn't find (recent enough) FFTWF. object not operational!", s->s_name);
#else
    pd_error(x, "[%s] compiled without FFTWF. object not operational!", s->s_name);
#endif
    warn_fftwf = 0;
  }


  x->x_setmultiout = (mtx_convolver_tilde_class == selected_class)?0:setmultiout;
  if (!x->x_setmultiout) {
    if (argc < 1) {
      x->ins = x->outs = 1;
    } else if (argc < 2) {
      x->ins = x->outs = (int)atom_getfloat(argv);
    } else {
      x->ins = (int)atom_getfloat(argv);
      x->outs = (int)atom_getfloat(argv + 1);
    }
    x->ins = (x->ins < 1)?1:x->ins;
    x->outs = (x->outs < 1)?1:x->outs;
    for (int i = 0; i < x->ins; i++) {
      inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    }
    for (int i = 0; i < x->outs; i++) {
      outlet_new(&x->x_obj, &s_signal);
    }
  }
  else {
    x->ins = x->outs = 0;
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
  x->coherent_xfade = 1;

  for (int i=1; i<argc; i++) { // use first symbol after argv[0] as init file name
    if (argv[i].a_type == A_SYMBOL) {
      t_symbol *matrix_file = argv[i].a_w.w_symbol;
      mtx_convolver_tilde_read(x, matrix_file);
      break;
    }
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

#ifdef HAVE_FFTWF
  my_functions.malloc = iemmatrix_get_stub("fftwf_malloc", mtx_convolver_tilde_class);
  my_functions.free = iemmatrix_get_stub("fftwf_free", mtx_convolver_tilde_class);
  my_functions.destroy_plan = iemmatrix_get_stub("fftwf_destroy_plan", mtx_convolver_tilde_class);
  my_functions.execute = iemmatrix_get_stub("fftwf_execute", mtx_convolver_tilde_class);
  my_functions.plan_dft_r2c_1d = iemmatrix_get_stub("fftwf_plan_dft_r2c_1d", mtx_convolver_tilde_class);
  my_functions.plan_dft_c2r_1d = iemmatrix_get_stub("fftwf_plan_dft_c2r_1d", mtx_convolver_tilde_class);
  IEMCONVOLVE(array_set_fftwf_functions) (&my_functions);

  warn_fftwf = !(IEMCONVOLVE(convolver_set_fftwf_functions) (&my_functions));
#endif
}
void iemtx_convolver__setup(void)
{
  mtx_convolver_tilde_setup();
}
