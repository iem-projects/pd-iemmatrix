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
  t_float x_pan;
  t_float f;
  int ins;
  int outs;
  int inlet_ins;
  t_float **hin;
  t_float **input;
  t_float **output;
  t_sample **inout_buffers;
  int rows;
  int cols;
  int m_inputs;
  int m_outputs;
  int m_irlen;
  t_canvas *x_canvas;
  conv_data *conv;
  t_inlet **x_in;
  t_outlet **x_out;
  int blocksize;
  int partition;
  int multichannel_mode;
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
  int min_ins = (x->m_inputs < x->ins) ? x->m_inputs : x->ins;
  int min_outs = (x->m_outputs < x->outs) ? x->m_outputs : x->outs;
  if (x->conv != 0) { // Operation: copy input signals, convolve, copy output signals
    for (int i = 0; i < min_ins; i++) {
      t_float *in = x->input[i];
      t_sample *pd_in = x->inout_buffers[i];
      for (int n = 0; n < x->blocksize; n++) {
        in[n] = (t_sample)pd_in[n];
      }
    }
    for (int i=min_ins; i < x->m_outputs; i++) { //necessary zeropad spurious channels case ??
      t_float *in = x->input[i];
      for (int n = 0; n < x->blocksize; n++) {
        in[n] = (t_sample)0;
      }
    }
    if (getNewIR(x->conv))
      post("perf2: have to crossfade to a new IR");
    conv_process(x->conv, x->input, x->output);
    for (int i = 0; i <min_outs; i++) {
      t_float *out = x->output[i];
      t_sample *pd_out = x->inout_buffers[i + x->ins];
      for (int n = 0; n < x->blocksize; n++) {
        pd_out[n] = (t_sample)out[n];
      }
    }
  } else { // DEFAULT: No convolver, deliver zero output signals.
    for (int i = 0; i < min_outs; i++) {
      t_sample *pd_out = x->inout_buffers[i + x->ins];
      for (int n = 0; n < x->blocksize; n++) {
        pd_out[n] = (t_sample)0;
      }
    }
  }
  return (w + 2);
}

void mtx_convolver_init(t_mtx_convolver_tilde *x, int ir_len) {
  if (x->conv)
    freeConvolution(x->conv);
  x->conv = 0;
  if ((x->blocksize > 0) && (x->ins > 0) && (x->outs > 0)) {
    int P = ceildiv(ir_len, x->blocksize);
    x->conv = initConvolution(x->blocksize, P, x->blocksize, x->ins, x->outs);
  } 
}

void mtx_convolver_resize(t_mtx_convolver_tilde *x, int new_blocksize) {
  int P = ceildiv(x->m_irlen,new_blocksize);
  if (x->conv) {
    if ((x->blocksize!=new_blocksize)||(x->conv->P!=P)||(x->conv->num_inputs!=x->m_inputs)
        ||(x->conv->num_outputs!=x->m_outputs)) {
      if (x->output)
        free2DArray(x->output, x->conv->num_outputs);
      if (x->input)
        free2DArray(x->input, x->conv->num_inputs);
      x->input = new2DArray(x->m_inputs, new_blocksize);
      x->output = new2DArray(x->m_outputs, new_blocksize);
      x->blocksize = new_blocksize;
    }
    freeConvolution(x->conv);
  }
  int ir_len = ceildiv(x->cols, x->blocksize) * x->blocksize;
  mtx_convolver_init(x, ir_len);
  setImpulseResponse2DZeropad(x->conv, x->hin, x->rows, x->cols);
  post("convolver/buffer new/resize: x->conv=%d, "
        "x->input=%d. x->output=%d",
        x->conv, x->input, x->output);
}

void mtx_convolver_tilde_dsp(t_mtx_convolver_tilde *x, t_signal **sp) {
  int ins = x->ins;
  int outs = x->outs;
  if (x->rows == 0 || x->cols == 0) {
    post("Warning: no valid matrix, skipping signal processing");
    // return;
  }
  if (x->multichannel_mode) { // override with num input signals in MC mode
    ins = sp[0]->s_nchans;
  }
  if (x->multichannel_mode) { // override with num input signals in MC mode
     // set number of outs in MC mode
    signal_setmultiout(&sp[1], x->m_outputs);
    outs = x->m_outputs;
  }
  x->ins = x->m_inputs;
  if ((ins+outs != x->ins+x->outs)){
    if (x->ins + x->outs > 0) {
      if (x->inout_buffers) {
        free(x->inout_buffers);
      }
      x->inout_buffers = (t_float **)malloc(sizeof(float *) * (ins + outs));
    }
  }
  x->outs = outs;
  x->ins = ins;
  int length = sp[0]->s_n;
  //mtx_convolver_resize(x,length);
  if (x->multichannel_mode) {
    for (int i = 0; i < x->ins; i++) {
      x->inout_buffers[i] = sp[0]->s_vec + i * length;
    }
    for (int i = 0; i < x->outs; i++) {
      x->inout_buffers[x->ins + i] = sp[1]->s_vec + i * length;
    }
  } else {
    for (int i = 0; i < x->ins + x->outs; i++) {
      x->inout_buffers[i] = sp[i]->s_vec;
    }
  }
  dsp_add(mtx_convolver_tilde_perform, 1, x);
}

void mtx_convolver_tilde_free(t_mtx_convolver_tilde *x) {
  if (x->inout_buffers)
    free(x->inout_buffers);
  if (x->hin)
    free2DArray(x->hin, x->rows);
  if (x->conv)
    freeConvolution(x->conv);
  if (x->input)
    free2DArray((t_float **)x->input, x->ins);
  if (x->output)
    free2DArray((t_float **)x->output, x->outs);
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
  x->input = 0;
  x->output = 0;
  x->inout_buffers = 0;
  x->hin = 0;
  x->conv = 0;
  x->blocksize = 0;
  x->rows = 0;
  x->cols = 0;

  return (void *)x;
}

const char*mtx_convolver_objname(void*obj) {
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
    pd_error(x, "%scrippled matrix", objname);
    return 1;
  }

  if ((tests & IEMMATRIX_CHECK_DIMENSIONS) && ((inputs<1)||(outputs<1)||(length<1))) {
    pd_error(x, "%sinvalid dimensions %dx%dx%d", objname, inputs, outputs, length);
    return 1;
  }

  if ((tests & IEMMATRIX_CHECK_SPARSE)&&(inputs* outputs*length>argc-3)) {
    pd_error(x, "%ssparse matrix not yet supported", objname);
    return 1;
  }
  return 0;
}

void mtx_convolver_tilde_matrix3D(t_mtx_convolver_tilde *x, t_symbol *s, int argc,
                      t_atom *argv) {
  int m_inputs, m_outputs, m_ir_len, ir_len = 0;
  
  if (argc < 3) {
    post("matrix3D message must have at least 3 arguments: inputs, outputs, ir_len");
    return;
  }
  if(mtx_convolver_check(x, argc, argv, 0))return;
  m_inputs = (int)atom_getfloat(argv);
  m_outputs = (int)atom_getfloat(argv + 1);
  m_ir_len = (int)atom_getfloat(argv + 2);
  argv += 3;
  if ((m_inputs == 0) || (m_outputs == 0) || (m_ir_len==0)) {
    post("%d inputs, %d outputs and %d ir_len specified", m_inputs, m_outputs, m_ir_len);
    return;
  } 

  post("inputs=%d, outputs=%d, ir_len=%d, x->blocksize=%d", m_inputs, m_outputs, m_ir_len, x->blocksize);

  int rows = m_inputs*m_outputs; int cols = m_ir_len; 
  if ((rows != x->rows) || (cols != x->cols)) {
    if (x->hin)
      free2DArray(x->hin, x->rows);
    x->hin = new2DArray(rows, cols);
    x->rows = rows;
    x->cols = cols;
    post("input mtx: re-sizing/setting x->rows=%d, x->cols=%d", x->rows,
         x->cols);
  }
  for (int io_idx = 0; io_idx < rows; io_idx++) { // store input matrix
    for (int time_index = 0; time_index < cols; time_index++) {
      x->hin[io_idx][time_index] = (t_float)atom_getfloat(argv++);
    }
  }

  if (x->blocksize) { // if blocksize is already known
      ir_len = ceildiv(cols, x->blocksize) * x->blocksize;
      mtx_convolver_init(x, ir_len);
    if (x->conv) { // if convolver exists: update IRs
      setImpulseResponse2DZeropad(x->conv, x->hin, x->rows, x->cols);
    }
  }
  x->m_inputs = m_inputs; 
  x->m_outputs= m_outputs;
}

static void mtx_convolver_tilde_read(t_mtx_convolver_tilde *x, t_symbol *filename)
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
      (!strcmp(ap->a_w.w_symbol->s_name,"matrix3D")
       || !strcmp(ap->a_w.w_symbol->s_name,"#matrix3D")) ) {
    mtx_convolver_tilde_matrix3D(x, gensym("matrix3D"), n, ap+1);
  }
  binbuf_free(bbuf);
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
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_matrix3D,
                  gensym("matrix3D"), A_GIMME, 0);
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_matrix3D,
  gensym("#matrix3D"), A_GIMME, 0);
  // class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_matrix,
  // gensym("matrix"), A_GIMME, 0);

  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_matrix3D,
                  gensym("matrix3D"), A_GIMME, 0);
  // class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_matrix,
  // gensym("matrix"), A_GIMME, 0);

  class_addmethod  (mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_read , gensym("read") ,
  A_SYMBOL, 0);
  class_addmethod  (mtx_convolver_tilde_class, (t_method)mtx_convolver_tilde_read , gensym("read") ,
  A_SYMBOL, 0);


  // CLASS_MAINSIGNALIN(mtx_convolver_tilde_class, t_mtx_convolver_tilde, f);
}
