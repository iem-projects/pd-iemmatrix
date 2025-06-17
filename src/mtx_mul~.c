/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) Thomas Musil, IEM KUG Graz Austria 2000-2003
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute, IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */


/* NOTES on multichannel
   - the number of incoming channels is not changeable by us
   - non-multichannel ports are always *hard* restrictions (wrong matrix is discarded)
   - the matrix dimension can only change while the DSP is off
     - while the DSP is on, non-matching matrices are discarded
   - if the dimensions change, the source and target matrices are resized
 */

#include "iemmatrix.h"

#ifndef CLASS_MULTICHANNEL
# define CLASS_MULTICHANNEL 0
#endif

typedef void (*setmultiout_f)(t_signal **sig, int nchans);


/* ---------- mtx_*~ - signal matrix multiplication object with message matrix-coeff. ----------- */


/* usage:
 *   [mtx_*~ <#out> <#in>]: multiply #in signals with a <#out,#in>-matrix to get #out-signals
 *   [mtx_*~ <#io>]: multiply #io signals with a <#io,#io>-matrix to get #io-signals
 *
 *  1st inlet: (message only): the matrix to multiply with (defaults to zero(#out, #in)
 *  2nd-(n-1)th inlet: the input signals
 *  last inlet: the interpolation time in [ms]
 *
 * this also implements a compatibility layer with old objects from tm and jmz
 * that provide basically the same functionality
 *
 * further note:
 *  pd does not like signal-objects that don't listen to signals on the very first inlet
 *  however, the first signal-inlet of [mtx_*~] is the 2nd(!) inlet
 *  to achieve this, we silently ignore any signal that comes in at the first inlet
 */

typedef struct _proxy {
  t_object p_obj;
  struct matrix_multilde*p_owner;
} t_proxy;

enum {
  IEMMATRIX = 0,
  IEM_MATRIX = 1, /* matrix_mul_line~ */
  ZEXY = 2        /* matrix~ */
} e_compat;


typedef struct matrix_multilde {
  /* private weirdo stuff at the beginning */
  t_object	x_obj;
  t_symbol      *x_name;
  t_proxy       *x_proxy;
  int           x_compat; /* 0=mtx_*~; 1=matrix_mul_line~; 2=matrix~ */
  setmultiout_f x_setmultiout; /* when doing multichannel, this is Pd>=0.54's signal_setmultiout(); otherwise NULL */

  /* DSP meta information */
  t_sample	**x_io; /* input/output signals (for easier passing to perform()) */
  t_float	x_msi; /* CLASS_MAINSIGNALIN() */
  int           x_dsp; /* is the DSP running? */

   /* number of non-multichannels iolets.
    * a value of '0' indicates a single multichannel iolets
    */
  size_t        x_inports, x_outports;

  /* the matrix */
  /* while DSP is running, rows=outs, cols=ins;
   * if XXports>0, then the rows/cols must match when setting a matrix
   */
  size_t        x_rows, x_cols;
  t_float       *x_matcur; /* current matrix (being interpolated) */
  t_float	*x_matend; /* interpolation target */

  /* interpolation data */
  t_float	x_time_ms; /* interpolation time (when new matrix comes through) */
  int		x_remaining_ticks; /* how long do we still need to interpolate */
  t_float	*x_inc, *x_biginc; /* matrices incrementation values */
  int		x_retarget; /* bool: whether we need to start a new interpolation */
  t_float	x_ms2tick; /* helper to translate time to ticks */

  /* buffer for matrix multiplication */
  t_sample	*x_outsumbuf; /* N samples for summing up */
  size_t	x_outsumbufsize;
} t_matrix_multilde;

static void proxy_dspstopped(t_proxy*p) {
  p->p_owner->x_dsp = 0;
}

static t_float*get_resized_matrix(const t_float*src,
				  const unsigned int srcR, const unsigned int srcC,
				  const unsigned int dstR, const unsigned int dstC) {
  t_float*dst = (t_float*)getbytes(dstR*dstC*sizeof(*dst));
  unsigned int c, r, cols, rows;
  if(!dst)
    return dst;
  cols = (srcC<dstC)?srcC:dstC;
  rows = (srcR<dstR)?srcR:dstR;
  for(r=0; r<rows; r++) {
    for(c=0; c<cols; c++) {
      dst[dstC*r + c] = src[srcC*r + c];
    }
  }
  return dst;
}
static t_float*resize_and_free(t_float*src,
			       const unsigned int srcR, const unsigned int srcC,
			       const unsigned int dstR, const unsigned int dstC) {
  t_float*dummy = get_resized_matrix(src, srcR, srcC, dstR, dstC);
  freebytes(src, srcR*srcC*sizeof(*src));
  return dummy;
}

static int matrix_multilde_resize(t_matrix_multilde *x, unsigned int rows, unsigned int cols) {
  if((rows == x->x_rows) && (cols == x->x_cols))
    return 1;
  x->x_matend = resize_and_free(x->x_matend, x->x_rows, x->x_cols, rows, cols);
  x->x_matcur = resize_and_free(x->x_matcur, x->x_rows, x->x_cols, rows, cols);
  x->x_inc    = resize_and_free(x->x_inc   , x->x_rows, x->x_cols, rows, cols);
  x->x_biginc = resize_and_free(x->x_biginc, x->x_rows, x->x_cols, rows, cols);
  if(!x->x_matend || !x->x_matcur || !x->x_inc || !x->x_biginc) {
    pd_error(x, "[%s] failed to resize matrices to [%dx%d]", x->x_name->s_name, rows, cols);
    x->x_rows = x->x_cols = 0;
    return 0;
  }
  x->x_rows = rows;
  x->x_cols = cols;
  return 1;
}


static t_class *matrix_multilde_mclass;
static t_class *matrix_multilde_class;
static t_class *matrix_multilde_proxy;

static void matrix_multilde_time(t_matrix_multilde *x, t_floatarg time_ms)
{
  if(time_ms <= 0.0f) {
    time_ms = 0.0f;
  }
  x->x_time_ms = time_ms;
}

static void matrix_multilde_matrix_set(t_matrix_multilde *x, int argc,
                                       t_atom *argv, int transpose)
{
  unsigned int col, row, length;
  t_float *matcur, *matend;

  if(argc<2) {
    pd_error(x, "[%s]: bad matrix: <int:out_rows> <int:in_cols> !", x->x_name->s_name);
    return;
  }
  if(iemmatrix_check(x, 0, argc, argv, 0))return;

  row = atom_getint(argv);
  argv++;
  col = atom_getint(argv);
  argv++;
  argc-=2;

  if(transpose) {
    int dummy=row;
    row=col;
    col=dummy;
  }

  if (x->x_dsp) {
    /* matrix dimensions must match while DSP is running */
    if ((col != x->x_cols) || (row != x->x_rows)) {
      pd_error(x, "[%s]: matrix dimensions must not change (%dx%d != %dx%d) while DSP is running!!",
	       x->x_name->s_name,
	       row, col, (int)x->x_rows, (int)x->x_cols);
      return;
    }
  } else {
    /* DSP is not running, check if we have a fixed number of iolets */
    if(x->x_inports && x->x_inports != col) {
      pd_error(x, "[%s]: cannot change fixed number of input channels (%d) to %d",
	       x->x_name->s_name, (int)x->x_inports, col);
      return;
    }
    if(x->x_outports && x->x_outports != row) {
      pd_error(x, "[%s]: cannot change fixed number of output channels (%d) to %d",
	       x->x_name->s_name, (int)x->x_outports, row);
      return;
    }

    if(!matrix_multilde_resize(x, row, col))
      return;
  }

  matcur = x->x_matcur;
  matend = x->x_matend;
  length = col * row;

  if(transpose) {
    /* we need to transpose the matrix */
    unsigned int i;
    for(i=0; i<row; i++) {
      unsigned int j;
      for(j=0; j<col; j++) {
        *matend++=atom_getfloat(argv+i+j*row);
      }
    }
  } else {
    unsigned int i;
    for(i=0; i<length; i++) {
      *matend++=atom_getfloat(argv++);
    }
  }

  if(x->x_time_ms <= 0.0f) {
    unsigned int i;
    matend = x->x_matend;
    for(i=0; i<length; i++) {
      *matcur++=*matend++;
    }
    x->x_remaining_ticks = x->x_retarget = 0;
  } else {
    x->x_retarget = 1;
  }
}
static void matrix_multilde_matrix(t_matrix_multilde *x, t_symbol *s,
                                   int argc, t_atom *argv)
{
  (void)s; /* unused */
  matrix_multilde_matrix_set(x, argc, argv, 0);
}
static void matrix_multilde_matrixT(t_matrix_multilde *x, t_symbol *s,
                                    int argc, t_atom *argv)
{
  /* transpose the matrix before setting it */
  (void)s; /* unused */
  matrix_multilde_matrix_set(x, argc, argv, 1);
}
static void matrix_multilde_element(t_matrix_multilde *x, t_symbol *s,
                                    int argc, t_atom *argv)
{
  unsigned int col, row, n_in_cols=x->x_cols;
  t_float element;
  t_float *matcur = x->x_matcur;
  t_float *matend = x->x_matend;
  (void)s; /* unused */

  if(argc != 3) {
    pd_error(x, "[%s]: bad arguments, expected <int:row> <int:column> <float:value>!", x->x_name->s_name);
    return;
  }

  int r = (int)atom_getfloat(argv+0) - 1;
  int c = (int)atom_getfloat(argv+1) - 1;
  row = (r >= 0)?r:0;
  col = (c >= 0)?c:0;

  element = atom_getfloat(argv+2);

  if((row >= x->x_rows) || (r < 0)) {
    pd_error(x, "[%s]: out of bound row!!", x->x_name->s_name);
    return;
  }
  if((col >= n_in_cols) || (c < 0)) {
    pd_error(x, "[%s]: out of bound column!!", x->x_name->s_name);
    return;
  }

  matend += row * n_in_cols + col;
  matcur += row * n_in_cols + col;

  if(x->x_time_ms <= 0.0f) {
    *matend = *matcur = element;
    x->x_remaining_ticks = x->x_retarget = 0;
  } else {
    *matend = element;
    x->x_retarget = 1;
  }
}

static void matrix_multilde_row(t_matrix_multilde *x, t_symbol *s,
                                int argc, t_atom *argv)
{
  unsigned int col, nth_row, i;
  t_float *matcur = x->x_matcur;
  t_float *matend = x->x_matend;
  (void)s; /* unused */

  if(argc<1) {
    pd_error(x,"[%s]: bad row!", x->x_name->s_name);
    return;
  }

  int _nth_row = (int)atom_getfloat(argv) -1;
  nth_row = (_nth_row>0)?_nth_row:0;
  argv++;
  argc--;

  if((nth_row >= x->x_rows) || (_nth_row < 0)) {
    pd_error(x, "[%s]: out of bound row!!", x->x_name->s_name);
    return;
  }
  col = x->x_cols;
  if(argc < 0 || ((unsigned int)argc < col)) {
    pd_error(x,"[%s]: col dimensions do not match !!", x->x_name->s_name);
    return;
  }

  matend += nth_row * col;
  matcur += nth_row * col;
  if(x->x_time_ms <= 0.0f) {
    for(i=0; i<col; i++) {
      *matend++ = *matcur++ = atom_getfloat(argv);
      argv++;
    }
    x->x_remaining_ticks = x->x_retarget = 0;
  } else {
    for(i=0; i<col; i++) {
      *matend++ = atom_getfloat(argv);
      argv++;
    }
    x->x_retarget = 1;
  }
}

static void matrix_multilde_col(t_matrix_multilde *x, t_symbol *s,
                                int argc, t_atom *argv)
{
  int row, col, nth_col, i;
  t_float *matcur = x->x_matcur;
  t_float *matend = x->x_matend;
  (void)s; /* unused */

  if(argc<1) {
    pd_error(x,"[%s]: bad column!", x->x_name->s_name);
    return;
  }

  nth_col = atom_getint(argv) - 1;
  argv++;
  argc--;

  col = x->x_cols;
  if((nth_col < 0)||(nth_col >= col)) {
    pd_error(x, "[%s]: out of bound column!!", x->x_name->s_name);
    return;
  }
  row = x->x_rows;
  if(argc < row) {
    pd_error(x,"[%s]: row dimensions do not match !!", x->x_name->s_name);
    return;
  }

  matend += nth_col;
  matcur += nth_col;
  if(x->x_time_ms <= 0.0f) {
    for(i=0; i<row; i++) {
      *matend = *matcur = atom_getfloat(argv);
      argv++;
      matend += col;
      matcur += col;
    }
    x->x_remaining_ticks = x->x_retarget = 0;
  } else {
    for(i=0; i<row; i++) {
      *matend = atom_getfloat(argv);
      argv++;
      matend += col;
      matcur += col;
    }
    x->x_retarget = 1;
  }
}

static void matrix_multilde_stop(t_matrix_multilde *x)
{
  int i = x->x_rows*x->x_cols;
  t_float *matend=x->x_matend;
  t_float *matcur=x->x_matcur;

  while(i--) {
    *matend++ = *matcur++;
  }
  x->x_remaining_ticks = x->x_retarget = 0;
}

/* the dsp thing */

static t_int *matrix_multilde_perf8(t_int *w)
{
  t_matrix_multilde *x = (t_matrix_multilde *)(w[1]);
  int n = (int)(w[2]);
  t_sample **io = x->x_io;
  t_sample *outsum, *houtsum;
  t_float *matcur, *matend;
  t_float *inc1 ,*biginc, inc;
  int n_in = x->x_cols;
  int n_out = x->x_rows;
  t_float *in, *out, mul, bigmul;
  int r, c, i;

  if(x->x_retarget) {
    int nticks = (int)(x->x_time_ms * x->x_ms2tick);

    if(!nticks) {
      nticks = 1;
    }
    x->x_remaining_ticks = nticks;
    inc1 = x->x_inc;
    biginc = x->x_biginc;
    matcur = x->x_matcur;
    matend = x->x_matend;
    bigmul = 1.0f / (t_float)nticks;
    mul = bigmul / ((t_float)n);
    i = n_out*n_in;
    while(i--) {
      inc = *matend++ - *matcur++;
      *inc1++ = inc * mul;
      *biginc++ = inc * bigmul;
    }
    x->x_retarget = 0;
  }

  if(x->x_remaining_ticks) {
    /* keep interpolating between current and target matrix */
    inc1 = x->x_inc;
    biginc = x->x_biginc;
    matcur = x->x_matcur;
    /* 1. output-vector-row */
    in = io[0];
    houtsum = x->x_outsumbuf;
    outsum = houtsum;
    inc = *inc1++;
    mul = *matcur;
    for(i=n; i; i -= 8, outsum += 8, in += 8) {
      outsum[0] = in[0] * mul;
      mul += inc;
      outsum[1] = in[1] * mul;
      mul += inc;
      outsum[2] = in[2] * mul;
      mul += inc;
      outsum[3] = in[3] * mul;
      mul += inc;
      outsum[4] = in[4] * mul;
      mul += inc;
      outsum[5] = in[5] * mul;
      mul += inc;
      outsum[6] = in[6] * mul;
      mul += inc;
      outsum[7] = in[7] * mul;
      mul += inc;
    }
    *matcur++ += *biginc++;
    for(c=1; c<n_in; c++) { /* c+1. element of 1. row */
      in = io[c];
      outsum = houtsum;
      inc = *inc1++;
      mul = *matcur;
      for(i=n; i; i -= 8, outsum += 8, in += 8) {
        outsum[0] += in[0] * mul;
        mul += inc;
        outsum[1] += in[1] * mul;
        mul += inc;
        outsum[2] += in[2] * mul;
        mul += inc;
        outsum[3] += in[3] * mul;
        mul += inc;
        outsum[4] += in[4] * mul;
        mul += inc;
        outsum[5] += in[5] * mul;
        mul += inc;
        outsum[6] += in[6] * mul;
        mul += inc;
        outsum[7] += in[7] * mul;
        mul += inc;
      }
      *matcur++ += *biginc++;
    }
    for(r=1; r<n_out; r++) { /* 2. .. n_out. output-vector-row */
      in = io[0];
      houtsum += n;
      outsum = houtsum;
      inc = *inc1++;
      mul = *matcur;
      for(i=n; i; i -= 8, outsum += 8, in += 8) {
        outsum[0] = in[0] * mul;
        mul += inc;
        outsum[1] = in[1] * mul;
        mul += inc;
        outsum[2] = in[2] * mul;
        mul += inc;
        outsum[3] = in[3] * mul;
        mul += inc;
        outsum[4] = in[4] * mul;
        mul += inc;
        outsum[5] = in[5] * mul;
        mul += inc;
        outsum[6] = in[6] * mul;
        mul += inc;
        outsum[7] = in[7] * mul;
        mul += inc;
      }
      *matcur++ += *biginc++;
      for(c=1; c<n_in; c++) { /* c+1. element of r+1. row */
        in = io[c];
        outsum = houtsum;
        inc = *inc1++;
        mul = *matcur;
        for(i=n; i; i -= 8, outsum += 8, in += 8) {
          outsum[0] += in[0] * mul;
          mul += inc;
          outsum[1] += in[1] * mul;
          mul += inc;
          outsum[2] += in[2] * mul;
          mul += inc;
          outsum[3] += in[3] * mul;
          mul += inc;
          outsum[4] += in[4] * mul;
          mul += inc;
          outsum[5] += in[5] * mul;
          mul += inc;
          outsum[6] += in[6] * mul;
          mul += inc;
          outsum[7] += in[7] * mul;
          mul += inc;
        }
        *matcur++ += *biginc++;
      }
    }

    if(!--x->x_remaining_ticks) {
      matcur = x->x_matcur;
      matend = x->x_matend;
      i = n_in * n_out;
      while(i--) {
        *matcur++ = *matend++;
      }
    }
  } else {
    /* no remaining ticks left - matrix stays constant for now */
    matend = x->x_matend;
    /* 1. output-vector-row */
    in = io[0];
    houtsum = x->x_outsumbuf;
    outsum = houtsum;
    mul = *matend++;
    if(mul == 0.0f) {
      for(i=n; i; i -= 8, outsum += 8, in += 8) {
        outsum[0] = 0.0f;
        outsum[1] = 0.0f;
        outsum[2] = 0.0f;
        outsum[3] = 0.0f;
        outsum[4] = 0.0f;
        outsum[5] = 0.0f;
        outsum[6] = 0.0f;
        outsum[7] = 0.0f;
      }
    } else {
      for(i=n; i; i -= 8, outsum += 8, in += 8) {
        outsum[0] = in[0] * mul;
        outsum[1] = in[1] * mul;
        outsum[2] = in[2] * mul;
        outsum[3] = in[3] * mul;
        outsum[4] = in[4] * mul;
        outsum[5] = in[5] * mul;
        outsum[6] = in[6] * mul;
        outsum[7] = in[7] * mul;
      }
    }
    for(c=1; c<n_in; c++) { /* c+1. element of 1. row */
      in = io[c];
      outsum = houtsum;
      mul = *matend++;
      if(mul != 0.0f) {
        for(i=n; i; i -= 8, outsum += 8, in += 8) {
          outsum[0] += in[0] * mul;
          outsum[1] += in[1] * mul;
          outsum[2] += in[2] * mul;
          outsum[3] += in[3] * mul;
          outsum[4] += in[4] * mul;
          outsum[5] += in[5] * mul;
          outsum[6] += in[6] * mul;
          outsum[7] += in[7] * mul;
        }
      }
    }
    for(r=1; r<n_out; r++) { /* 2. .. n_out. output-vector-row */
      in = io[0];
      houtsum += n;
      outsum = houtsum;
      mul = *matend++;
      if(mul == 0.0f) {
        for(i=n; i; i -= 8, outsum += 8, in += 8) {
          outsum[0] = 0.0f;
          outsum[1] = 0.0f;
          outsum[2] = 0.0f;
          outsum[3] = 0.0f;
          outsum[4] = 0.0f;
          outsum[5] = 0.0f;
          outsum[6] = 0.0f;
          outsum[7] = 0.0f;
        }
      } else {
        for(i=n; i; i -= 8, outsum += 8, in += 8) {
          outsum[0] = in[0] * mul;
          outsum[1] = in[1] * mul;
          outsum[2] = in[2] * mul;
          outsum[3] = in[3] * mul;
          outsum[4] = in[4] * mul;
          outsum[5] = in[5] * mul;
          outsum[6] = in[6] * mul;
          outsum[7] = in[7] * mul;
        }
      }
      for(c=1; c<n_in; c++) { /* c+1. element of r+1. row */
        in = io[c];
        outsum = houtsum;
        mul = *matend++;
        if(mul != 0.0f) {
          for(i=n; i; i -= 8, outsum += 8, in += 8) {
            outsum[0] += in[0] * mul;
            outsum[1] += in[1] * mul;
            outsum[2] += in[2] * mul;
            outsum[3] += in[3] * mul;
            outsum[4] += in[4] * mul;
            outsum[5] += in[5] * mul;
            outsum[6] += in[6] * mul;
            outsum[7] += in[7] * mul;
          }
        }
      }
    }
  }
  outsum = x->x_outsumbuf;
  for(r=0; r<n_out; r++) { /* output-vector-row */
    out = io[n_in+r];
    for (i=n; i; i -= 8, out += 8, outsum += 8) {
      out[0] = outsum[0];
      out[1] = outsum[1];
      out[2] = outsum[2];
      out[3] = outsum[3];
      out[4] = outsum[4];
      out[5] = outsum[5];
      out[6] = outsum[6];
      out[7] = outsum[7];
    }
  }
  return (w+3);
}
static t_int *matrix_multilde_perform(t_int *w)
{
  t_matrix_multilde *x = (t_matrix_multilde *)(w[1]);
  int n = (int)(w[2]);
  t_sample **io = x->x_io;
  t_sample *outsum, *houtsum;
  t_float *matcur, *matend;
  t_float *inc1 ,*biginc, inc;
  int n_in = x->x_cols;
  int n_out = x->x_rows;
  t_float *in, *out, mul, bigmul;
  int r, c, i;

  if(x->x_retarget) {
    int nticks = (int)(x->x_time_ms * x->x_ms2tick);

    if(!nticks) {
      nticks = 1;
    }
    x->x_remaining_ticks = nticks;
    inc1 = x->x_inc;
    biginc = x->x_biginc;
    matcur = x->x_matcur;
    matend = x->x_matend;
    bigmul = 1.0f / ((t_float)nticks);
    mul = bigmul / ((t_float)n);
    i = n_out*n_in;
    while(i--) {
      inc = *matend++ - *matcur++;
      *inc1++ = inc * mul;
      *biginc++ = inc * bigmul;
    }
    x->x_retarget = 0;
  }

  if(x->x_remaining_ticks) {
    /* keep interpolating between current and target matrix */
    inc1 = x->x_inc;
    biginc = x->x_biginc;
    matcur = x->x_matcur;
    /* 1. output-vector-row */
    in = io[0];
    houtsum = x->x_outsumbuf;
    outsum = houtsum;
    inc = *inc1++;
    mul = *matcur;
    i=n;
    while(i--) {
      *outsum++ = *in++ * mul;
      mul += inc;
    }
    *matcur++ += *biginc++;
    for(c=1; c<n_in; c++) { /* c+1. element of 1. row */
      in = io[c];
      outsum = houtsum;
      inc = *inc1++;
      mul = *matcur;
      i=n;
      while(i--) {
        *outsum++ += *in++ * mul;
        mul += inc;
      }
      *matcur++ += *biginc++;
    }
    for(r=1; r<n_out; r++) { /* 2. .. n_out. output-vector-row */
      in = io[0];
      houtsum += n;
      outsum = houtsum;
      inc = *inc1++;
      mul = *matcur;
      i=n;
      while(i--) {
        *outsum++ = *in++ * mul;
        mul += inc;
      }
      *matcur++ += *biginc++;
      for(c=1; c<n_in; c++) { /* c+1. element of r+1. row */
        in = io[c];
        outsum = houtsum;
        inc = *inc1++;
        mul = *matcur;
        i=n;
        while(i--) {
          *outsum++ += *in++ * mul;
          mul += inc;
        }
        *matcur++ += *biginc++;
      }
    }

    if(!--x->x_remaining_ticks) {
      matcur = x->x_matcur;
      matend = x->x_matend;
      i = n_in * n_out;
      while(i--) {
        *matcur++ = *matend++;
      }
    }
  } else {
    /* no remaining ticks left - matrix stays constant for now */
    matend = x->x_matend;
    /* 1. output-vector-row */
    in = io[0];
    houtsum = x->x_outsumbuf;
    outsum = houtsum;
    mul = *matend++;
    i=n;
    if(mul == 0.0f)
      while(i--) {
        *outsum++ = 0.0f;
      }
    else
      while(i--) {
        *outsum++ = *in++ * mul;
      }

    for(c=1; c<n_in; c++) { /* c+1. element of 1. row */
      in = io[c];
      outsum = houtsum;
      mul = *matend++;
      if(mul != 0.0f) {
        i=n;
        while(i--) {
          *outsum++ += *in++ * mul;
        }
      }
    }
    for(r=1; r<n_out; r++) { /* 2. .. n_out. output-vector-row */
      in = io[0];
      houtsum += n;
      outsum = houtsum;
      mul = *matend++;
      i=n;
      if(mul == 0.0f)
        while(i--) {
          *outsum++ = 0.0f;
        }
      else
        while(i--) {
          *outsum++ = *in++ * mul;
        }

      for(c=1; c<n_in; c++) { /* c+1. element of r+1. row */
        in = io[c];
        outsum = houtsum;
        mul = *matend++;
        i=n;
        if(mul != 0.0f)
          while(i--) {
            *outsum++ += *in++ * mul;
          }
      }
    }
  }
  outsum = x->x_outsumbuf;

  for(r=0; r<n_out; r++) { /* output-vector-row */
    out = io[n_in+r];
    i=n;
    while(i--) {
      *out++ = *outsum++;
    }
  }

  return (w+3);
}


static void matrix_multilde_dsp(t_matrix_multilde *x, t_signal **sp)
{
  const unsigned int length = (sp[0]->s_n > 0)?sp[0]->s_n:0;
  unsigned int i, n=length * x->x_rows;
  /* [mtx_*~] ignores the signal on the very 1st inlet */
  int compat_offset=(x->x_compat)?0:1;
  size_t ichannels = x->x_inports?x->x_inports:0;
  size_t ochannels = x->x_outports?x->x_outports:x->x_rows;
  t_sample**io;

  /* DSP is running */
  x->x_dsp = 1;

  if(x->x_setmultiout) {
#if CLASS_MULTICHANNEL
    size_t inports = (x->x_inports < 1)?1:x->x_inports;
    /* multichannel mode */
    if(!ichannels) {
      /* multichannel input */
      ichannels = sp[compat_offset]->s_nchans;
    }

    if(x->x_outports) {
      /* create singlechannel outputs */
      for(i=0; i<ochannels; i++) {
	x->x_setmultiout(&sp[compat_offset + inports + i], 1);
      }
    } else {
      /* create multichannel output */
      ochannels = x->x_rows;
      x->x_setmultiout(&sp[compat_offset + inports], ochannels);
    }
#else
    pd_error(x, "BUG: multichannel enabled but not compile-time multichannel support");
    return;
#endif
  }

  /* ensure that the matrices have the correct dimension */
  io = (t_sample **)resizebytes(x->x_io,
      (x->x_cols + x->x_rows) * sizeof(*x->x_io),
      (ichannels + ochannels) * sizeof(*x->x_io));

  if(!io) {
    pd_error(x, "Unable to get memory");
    return;
  }
  x->x_io = io;

  if(!matrix_multilde_resize(x, ochannels, ichannels)) {
    return;
  }

  if(!x->x_outsumbuf) {
    x->x_outsumbufsize = n;
    x->x_outsumbuf = (t_sample*)getbytes(x->x_outsumbufsize * sizeof(t_sample));
  } else if(x->x_outsumbufsize != n) {
    x->x_outsumbuf = (t_sample*)resizebytes(x->x_outsumbuf,
                                            x->x_outsumbufsize*sizeof(t_sample),
                                            n*sizeof(t_sample));
    x->x_outsumbufsize = n;
  }

  if(x->x_setmultiout) {
#if CLASS_MULTICHANNEL
    int offset = compat_offset;
    /* setup input channels */
    if(x->x_inports) {
      for(i=0; i<ichannels; i++) {
	x->x_io[i] = sp[offset + i]->s_vec;
      }
      offset += ichannels;
    } else {
      for(i=0; i<ichannels; i++) {
	x->x_io[i] = sp[offset]->s_vec + i * length;
      }
      offset += 1;
    }
    /* setup output channels */
    if(x->x_outports) {
      for(i=0; i<ochannels; i++) {
	x->x_io[ichannels + i] = sp[offset + i]->s_vec;
      }
      offset += ichannels;
    } else {
      for(i=0; i<ochannels; i++) {
	x->x_io[ichannels + i] = sp[offset]->s_vec + i * length;
      }
      offset += 1;
    }
#endif
  } else {
    n = ichannels + ochannels;
    for(i=0; i<n; i++) {
      x->x_io[i] = sp[compat_offset+i]->s_vec;
    }
  }

  x->x_ms2tick = 0.001f * (t_float)(sp[0]->s_sr) / (t_float)length;


  if(length&7) {
    dsp_add(matrix_multilde_perform, 2, x, length);
  } else {
    dsp_add(matrix_multilde_perf8, 2, x, length);
  }
}


/* setup/setdown things */

static void matrix_multilde_free(t_matrix_multilde *x)
{
  freebytes(x->x_matcur, x->x_cols * x->x_rows * sizeof(t_float));
  freebytes(x->x_matend, x->x_cols * x->x_rows * sizeof(t_float));
  freebytes(x->x_inc, x->x_cols * x->x_rows * sizeof(t_float));
  freebytes(x->x_biginc, x->x_cols * x->x_rows * sizeof(t_float));
  freebytes(x->x_io, (x->x_cols + x->x_rows) * sizeof(t_sample*));
  if(x->x_outsumbuf) {
    freebytes(x->x_outsumbuf, x->x_outsumbufsize * sizeof(t_sample));
  }
  if(x->x_proxy) {
    pd_unbind(&x->x_proxy->p_obj.ob_pd, gensym("pd-dsp-stopped"));
    pd_free(&x->x_proxy->p_obj.ob_pd);
  }
}

static void *matrix_multilde_new(t_symbol *s, int argc, t_atom *argv)
{
  /*
    // singlechannel
    [mtx_*~       ] singlechannel in:1/1 out:1/1
    [mtx_*~ 3     ] singlechannel in:3/3 out:3/3
    [mtx_*~ 3 4   ] singlechannel in:4/4 out:3/3

    // in/out is singlechannel, out/in is multichannel
    [mtx_*~ 0 4   ] mixedchannel  in:4/4 out:1/?
    [mtx_*~ 3 0   ] mixedchannel  in:1/? out:3/3

    // implicit full multichannel
    [mtx_*~ 0     ] multichannel in:1/? out:1/?
    [mtx_*~ 0 0   ] multichannel in:1/? out:1/?
    // explicit full multichannel
    [mtx_*~ -m    ] multichannel in:1/? out:1/?
    [mtx_*~ -m 0  ] multichannel in:1/? out:1/?
    [mtx_*~ -m 0 0] multichannel in:1/? out:1/?
    [mtx_*~ -m 3  ] multichannel in:1/3 out:1/3
    [mtx_*~ -m 3 4] multichannel in:1/4 out:1/3 (useless, as inchannels are defined by signal-in)
    [mtx_*~ -m 0 4] multichannel in:1/4 out:1/? (useless, as inchannels are defined by signal-in)
    [mtx_*~ -m 3 0] multichannel in:1/? out:1/3 (useless, as inchannels are defined by signal-in)


    // other names
    [matrix_mul~] ALIAS
    [mtx_mul~]    ALIAS

    [matrix_mul_line~]  LEGACY (iem_matrix); always singlechannel
    [matrix~] LEGACY (zexy); always singlechannel
   */
  t_matrix_multilde *x = 0;
  t_class*cls;
  int force_multi = 0, want_multi = 0;
  setmultiout_f setmultiout = iemmatrix_getpdfun("signal_setmultiout");
  unsigned int portsIn, portsOut, channelsIn, channelsOut;
  int nin, nout;
  t_float interpoltime = 0.;
  int compat = IEMMATRIX;
  int i;

  t_atom*ap_in, *ap_out, *ap_time;

  if(s==gensym("matrix~")) {
    compat = ZEXY;
  } else if (s==gensym("matrix_mul_line~")) {
    compat = IEM_MATRIX;
  } else {
    if(argc && A_SYMBOL == argv->a_type) {
      if (atom_getsymbol(argv) == gensym("-m")) {
	force_multi = want_multi = 1;
	argc--;
	argv++;
      } else {
	goto usage;
      }
    }
  }
  if(argc && A_FLOAT != argv->a_type) {
    goto usage;
  }

  /* arguments parsing:
   *  this might depend on whether we are creating an object
   *  [mtx_*~], [matrix~] or [matrix_mul_line~]
   *
   * [mtx_*~  [#out [#in [time]]]]:  in1=matrix; in2-in(#in+1)=signals; in(#in+2)=time
   * [matrix~ [#in [#out [time]]]]:  in1-in(#in)=signals; in(#in+1)=matrix; in(#in+2)=time
   * [matrix_mul_line~ [#in [#out [time]]]]: in1=matrix; in1=time; in1-in(#in):=signals
   *
   * furthermore:
   *  [mtx_*~] and [matrix_mul_line~] : O^=A*I^
   *  [matrix~]                       : O^'=I^'*B
   *
   *  with "matrix=(A or B)" and "A=B'"
   */
  ap_in  =argv+(compat?0:1);
  ap_out =argv+(compat?1:0);
  ap_time=argv+2;
  switch(argc) {
  case 0:
    nin = nout = 1;
    interpoltime = (ZEXY == compat)?0.f:50.0f;
    break;
  case 1:
    nin = nout = (int)atom_getfloat(argv);
    interpoltime = (ZEXY == compat)?0.f:50.0f;
    break;
  case 2:
    nin = (int)atom_getfloat(ap_in);
    nout = (int)atom_getfloat(ap_out);
    interpoltime = (ZEXY == compat)?0.f:50.0f;
    break;
  default:
    nin = (int)atom_getfloat(ap_in);
    nout = (int)atom_getfloat(ap_out);
    interpoltime = atom_getfloat(ap_time);
    break;
  }

  /* make sure that these are unsigned */
  if((nin < 0) || (nout < 0)) {
    pd_error(0, "[%s] matrix dimensions must not be negative [%dx%d]", s->s_name, nout, nin);
    return 0;
  }

  /* sanity check */
  if(interpoltime < 0.0f) {
    interpoltime = (IEM_MATRIX == compat)?50.f:0.0f;
  }


  /* do we need multichannel? */
  if (!compat && (0==nin || 0==nout))
    want_multi = 1;

  /* can we do multichannel? */
  cls = matrix_multilde_class;
  if(!compat && want_multi) {
    if (CLASS_MULTICHANNEL && setmultiout)
      cls = matrix_multilde_mclass;
  }

  /* relate ports/channels to the arguments */
  portsIn = channelsIn = nin;
  portsOut = channelsOut = nout;
  if(force_multi) {
    portsIn = portsOut = 0;
  }
  if(channelsIn <1) channelsIn  = 1;
  if(channelsOut<1) channelsOut = 1;

  /* create the object */
  x = (t_matrix_multilde *)pd_new(cls);
  if(compat) {
    /* print an error, nudging the user to upgrade to [mtx_*~] */

    /* however, suppress the error in help-patches */
    const char*absname = iemmatrix_parentabstractionname(0);
    if(absname) {
      char*absext = strrchr(absname, '-');
      if (absext && !strcmp(absext, "-help.pd")) {
        /* the object lives in a help-patch: suppress error */
      } else {
        absname = 0;
      }
    }
    if(!absname)
      pd_error(x, "[%s] is deprecated! use [mtx_*~] instead!", s->s_name);
  }
  x->x_name = s;
  x->x_compat=compat;
  x->x_setmultiout = (matrix_multilde_class == cls)?0:setmultiout;

  if (!compat && want_multi) {
    /* user requested multichannel */
    /* warn (but only once) if we cannot actually do that */
    static int warn_multichannel = 1;
    if(!x->x_setmultiout) {
      int major, minor, bugfix;
      sys_getversion(&major, &minor, &bugfix);
      if(warn_multichannel)
        pd_error(x, "[%s] multichannel requested, but iemmatrix is running in Pd-%d.%d-%d, which doesn't support it", s->s_name, major, minor, bugfix);
      x->x_setmultiout = 0;
    } else if (!CLASS_MULTICHANNEL) {
      if(warn_multichannel)
        pd_error(x, "[%s] multichannel requested, but iemmatrix was built against Pd-%d.%d-%d, which doesn't support it", s->s_name, PD_MAJOR_VERSION, PD_MINOR_VERSION, PD_BUGFIX_VERSION);
      x->x_setmultiout = 0;
    }
    warn_multichannel = 0;
  }

  x->x_proxy = (t_proxy*)pd_new(matrix_multilde_proxy);
  x->x_proxy->p_owner = x;
  pd_bind(&x->x_proxy->p_obj.ob_pd, gensym("pd-dsp-stopped"));

  x->x_inports = portsIn;
  x->x_outports = portsOut;
  x->x_time_ms = interpoltime;

  x->x_rows = channelsOut;
  x->x_cols = channelsIn;

  /* creating signal ins & outs */
  /* in compat mode, the 1st signal inlet is already made */
  if(portsIn<1)portsIn=1;
  if(portsOut<1)portsOut=1;
  i = portsIn - !!x->x_compat;
  while(i--) {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  }
  i = portsOut;
  while(i--) {
    outlet_new(&x->x_obj, &s_signal);
  }

  /* creating the matrix-inlet for [matrix~] */
  if(ZEXY == x->x_compat) {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"), gensym(""));
  }

  /* creating time-inlet (not for [matrix_mul_linie~]) */
  if(x->x_compat != IEM_MATRIX) {
    inlet_new(&x->x_obj,  &x->x_obj.ob_pd, &s_float, gensym("time"));
  }

  /* setting up internal values */
  x->x_matcur = (t_float *)getbytes(x->x_cols * x->x_rows * sizeof(*x->x_matcur));
  x->x_matend = (t_float *)getbytes(x->x_cols * x->x_rows * sizeof(*x->x_matend));
  x->x_inc = (t_float *)getbytes(x->x_cols * x->x_rows * sizeof(*x->x_inc));
  x->x_biginc = (t_float *)getbytes(x->x_cols * x->x_rows * sizeof(*x->x_biginc));
  x->x_io = (t_sample **)getbytes((x->x_cols + x->x_rows) * sizeof(*x->x_io));

  return (x);

usage:
  pd_error(0, "[%s] bad arguments\n\tuse '<int:ouchannels> <int:inchannels>'\n\tor '-m <int:outchannels> <int:inchannels>'", s->s_name);
  return 0;
}

static void mtx_mul_addmethods(t_class*c) {
  class_addmethod(c, (t_method)matrix_multilde_dsp,
                  gensym("dsp"), 0);

  class_addmethod(c, (t_method)matrix_multilde_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addmethod(c, (t_method)matrix_multilde_element,
                  gensym("element"), A_GIMME, 0);
  class_addmethod(c, (t_method)matrix_multilde_row,
                  gensym("row"), A_GIMME, 0);
  class_addmethod(c, (t_method)matrix_multilde_col,
                  gensym("col"), A_GIMME, 0);
  class_addmethod(c, (t_method)matrix_multilde_stop,
                  gensym("stop"), 0);
  class_addmethod(c, (t_method)matrix_multilde_time,
                  gensym("time"), A_FLOAT, 0);

  /* LATER: can we reuse the 'matrix' method and transpose
   * depending on the compat level? */
  class_addmethod(c, (t_method)matrix_multilde_matrixT,
                  gensym(""), A_GIMME, 0);

  CLASS_MAINSIGNALIN(c, t_matrix_multilde, x_msi);
}

void mtx_mul_tilde_setup(void)
{
  if (CLASS_MULTICHANNEL && iemmatrix_getpdfun("signal_setmultiout")) {
    /* multichannel variant */
    matrix_multilde_mclass = class_new(gensym("mtx_mul~"),
				      (t_newmethod)matrix_multilde_new,
				      (t_method)matrix_multilde_free,
				      sizeof(t_matrix_multilde),
				      0 | CLASS_MULTICHANNEL,
				      A_GIMME, 0);
    /* non-multichannel variant */
    /* compatibility with jmz's zexy */
    matrix_multilde_class = class_new(gensym("matrix~"),
					     (t_newmethod)matrix_multilde_new,
					     (t_method)matrix_multilde_free,
					     sizeof(t_matrix_multilde),
					     0,
					     A_GIMME, 0);
    class_sethelpsymbol(matrix_multilde_class, gensym("mtx_mul~"));
  } else {
    matrix_multilde_class = class_new(gensym("mtx_mul~"),
				      (t_newmethod)matrix_multilde_new,
				      (t_method)matrix_multilde_free,
				      sizeof(t_matrix_multilde),
				      0,
				      A_GIMME, 0);
    class_addcreator((t_newmethod)matrix_multilde_new, gensym("matrix~"),
                     A_GIMME, 0);
    matrix_multilde_mclass = matrix_multilde_class;
  }

  class_addcreator((t_newmethod)matrix_multilde_new, gensym("mtx_*~"),
                   A_GIMME, 0);
  class_addcreator((t_newmethod)matrix_multilde_new, gensym("matrix_mul~"),
                   A_GIMME, 0);

  /* compatibility with tm's iem_matrix */
  class_addcreator((t_newmethod)matrix_multilde_new,
                   gensym("matrix_mul_line~"), A_GIMME, 0);

  mtx_mul_addmethods(matrix_multilde_mclass);
  if(matrix_multilde_mclass != matrix_multilde_class)
    mtx_mul_addmethods(matrix_multilde_class);

  matrix_multilde_proxy = class_new(gensym("mtx_*~ proxy"),
				    0, 0,
				    sizeof(t_proxy),
				    CLASS_PD,
				    0);
  class_addbang(matrix_multilde_proxy, (t_method)proxy_dspstopped);
}

void iemtx_mul__setup(void)
{
  mtx_mul_tilde_setup();
}
