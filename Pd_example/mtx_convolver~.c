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

static t_class *mtx_convolver_tilde_class;
static t_class *mtx_convolver_tilde_mclass;

typedef struct _mtx_convolver_tilde {
  t_object x_obj;
  t_float x_pan;
  t_float f;
  int ins;
  int outs;
  int renew_convolver;
  float **hin;
  float **input;
  float **output;
  t_sample **inout_buffers;
  int rows;
  int cols;
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
  if (x->conv !=
      0) { // Operation: copy input signals, convolve, copy output signals
    for (int i = 0; i < x->ins; i++) {
      float *in = x->input[i];
      t_sample *pd_in = x->inout_buffers[i];
      for (int n = 0; n < x->blocksize; n++) {
        in[n] = (float)pd_in[n];
      }
    }
    if (getNewIR(x->conv))
      post("perf2: have to crossfade to a new IR");
    conv_process(x->conv, x->input, x->output);
    for (int i = 0; i < x->outs; i++) {
      float *out = x->output[i];
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

void mtx_convolver_init(t_mtx_convolver_tilde *x, int ir_len) {
  if (x->conv)
    freeConvolution(x->conv);
  x->conv = 0;
  if ((x->blocksize > 0) && (x->ins > 0) && (x->outs > 0)) {
    int P = ceildiv(ir_len, x->blocksize);
    x->conv = initConvolution(x->blocksize, P, x->blocksize, x->ins, x->outs);
  }
}

void mtx_convolver_tilde_dsp(t_mtx_convolver_tilde *x, t_signal **sp) {
  int ins = x->ins;
  int outs = x->outs;
  if (x->rows == 0 || x->cols == 0) {
    post("Warning: no valid matrix, skipping signal processing");
  }
  if (x->multichannel_mode) { // override with num input signals in MC mode
    ins = sp[0]->s_nchans;
  }
  outs = ceildiv(x->rows, ins); // calculate number of outs from stored matrix
  if (x->multichannel_mode) {   // set number of outs in MC mode
    signal_setmultiout(&sp[1], outs);
  }
  if ((outs != x->outs) || (ins != x->ins)) {
    post("mtx_convolver~ DSP: resizing from %d ins %d outs to %d ins %d outs",
         x->ins, x->outs, ins, outs);
    if (x->output)
      free2DArray(x->output, x->outs);
    if (x->input)
      free2DArray(x->input, x->ins);
    if (x->inout_buffers)
      free(x->inout_buffers);
    x->input = 0;
    x->output = 0;
    x->outs = outs;
    x->ins = ins;
    x->inout_buffers = 0;
    if (x->ins > 0) {
      x->input = new2DArray(x->ins, sp[0]->s_n);
    }
    if (x->outs > 0) {
      x->output = new2DArray(x->outs, sp[0]->s_n);
    }
    if (x->ins + x->outs > 0) {
      x->inout_buffers = (float **)malloc(sizeof(float *) * (x->ins + x->outs));
    }
  }
  int length = sp[0]->s_n;
  x->renew_convolver = 1;
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
  if (x->blocksize != sp[0]->s_n) {
    if (x->output)
      free2DArray(x->output, x->outs);
    if (x->input)
      free2DArray(x->input, x->ins);
    x->input = new2DArray(x->ins, sp[0]->s_n);
    x->output = new2DArray(x->outs, sp[0]->s_n);
    x->blocksize = sp[0]->s_n;
    x->renew_convolver = 1;
  }
  if (x->renew_convolver) {
    int ir_len = ceildiv(x->cols, x->blocksize) * x->blocksize;
    mtx_convolver_init(x, ir_len);
    if (x->conv) {
      setImpulseResponse2DZeropad(x->conv, x->hin, x->rows, x->cols);
      post("dsp2: convolver did not exist previously\n new: x->conv=%d, "
           "x->input=%d. x->output=%d",
           x->conv, x->input, x->output);
    }
    x->renew_convolver = 0;
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
    free2DArray((float **)x->input, x->ins);
  if (x->output)
    free2DArray((float **)x->output, x->outs);
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
  x->ins = x->outs = 1;
  if (!x->multichannel_mode) {
    if (argc < 1) {
      x->ins = x->outs = 1;
    } else if (argc < 2) {
      x->ins = x->outs = (int)atom_getfloat(argv);
    } else {
      x->ins = (int)atom_getfloat(argv);
      x->outs = (int)atom_getfloat(argv + 1);
    }
  }

  if (x->multichannel_mode) {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
  } else {
    for (int i = 0; i < x->ins; i++) {
      inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    }
    for (int i = 0; i < x->outs; i++) {
      outlet_new(&x->x_obj, &s_signal);
    }
  }
  x->input = 0;
  x->output = 0;
  x->inout_buffers = 0;
  x->hin = 0;
  x->conv = 0;
  x->blocksize = 0;
  x->rows = 0;
  x->cols = 0;
  x->renew_convolver = 1;

  return (void *)x;
}

void mtx_input_matrix(t_mtx_convolver_tilde *x, t_symbol *s, int argc,
                      t_atom *argv) {
  int rows, cols, ir_len = 0;
  int mc_changed_outs = 0;
  if (argc < 3) {
    post("no valid matrix with just 2 entries");
    return;
  }
  rows = (int)atom_getfloat(argv);
  cols = (int)atom_getfloat(argv + 1);
  if (rows * cols > argc - 2) {
    post("no valid matrix: only %d elements instead of r x c=%d", argc - 2,
         rows * cols);
    return;
  }
  if ((cols == 0) || (rows == 0)) {
    post("%d columns, %d rows specified", cols, rows);
    return;
  }
  argv += 2;
  post("rows=%d, cols=%d, x->blocksize=%d", rows, cols, x->blocksize);
  if ((rows != x->rows) || (cols != x->cols)) {
    x->renew_convolver = 1;
    if (x->hin)
      free2DArray(x->hin, x->rows);
    x->hin = new2DArray(rows, cols);
    x->rows = rows;
    x->cols = cols;
    post("input mtx: re-sizing/setting x->rows=%d, x->cols=%d", x->rows,
         x->cols);
    if (x->multichannel_mode) {
      mc_changed_outs = x->outs;
      x->outs = ceildiv(x->rows, x->ins); // number of output signals
    }
  }
  for (int io_idx = 0; io_idx < rows; io_idx++) { // store input matrix
    for (int time_index = 0; time_index < cols; time_index++) {
      x->hin[io_idx][time_index] = (float)atom_getfloat(argv++);
    }
  }
  if (x->blocksize) { // if blocsize is already known
    if (x->renew_convolver) {
      ir_len = ceildiv(cols, x->blocksize) * x->blocksize;
      mtx_convolver_init(x, ir_len);
      if ((mc_changed_outs != x->outs) &&
          (canvas_dspstate)) { // outputs changed while dsp running
        post("DSP update to get %d outputs instead of %d", x->outs, mc_changed_outs);
        canvas_update_dsp();
      }
      x->renew_convolver = 0;
    }
    if (x->conv) { // if convolver exists: update IRs
      setImpulseResponse2DZeropad(x->conv, x->hin, x->rows, x->cols);
    }
  }
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
  class_addmethod(mtx_convolver_tilde_class, (t_method)mtx_input_matrix,
                  gensym("matrix"), A_GIMME, 0);

  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_convolver_tilde_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(mtx_convolver_tilde_mclass, (t_method)mtx_input_matrix,
                  gensym("matrix"), A_GIMME, 0);

  // CLASS_MAINSIGNALIN(mtx_convolver_tilde_class, t_mtx_convolver_tilde, f);
}