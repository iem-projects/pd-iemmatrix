#include "iemmatrix.h"
#define MTX_PACK_MAXCHANNELS 200

static t_class *mtx_unpack_tilde_class;

typedef struct _mtx_unpack_tilde {
  t_object x_obj;
  int rows;
  int cols;
  int block_size;
  int num_chan;
  t_float **sig_out;
  t_atom *list_in;
  t_int *(*perform_fcn)(t_int*);
} mtx_unpack_tilde;

static t_int *mtx_unpack_Perform (t_int *arg)
{
  mtx_unpack_tilde *x = (mtx_unpack_tilde *) (arg[1]);
  return (x->perform_fcn(arg));
}

static t_int *mtx_unpack_PerformInactive (t_int *arg)
{
  return(arg+2);
}

static t_int *mtx_unpack_PerformSetInactive (t_int *arg)
{
  mtx_unpack_tilde *x = (mtx_unpack_tilde *) (arg[1]);
  int chan;
  int samp;
  t_atom *lptr=x->list_in;

  for (chan=0; chan<x->num_chan; chan++) {
    for (samp=0; samp<x->block_size; samp++) {
      x->sig_out[chan][samp]=0;
    }
    lptr+=x->cols;
  }
  x->perform_fcn=mtx_unpack_PerformInactive;
  return(arg+2);
}

static t_int *mtx_unpack_PerformActive (t_int *arg)
{
  mtx_unpack_tilde *x = (mtx_unpack_tilde *) (arg[1]);
  int chan;
  int samp;
  const int maxchan = (x->rows < x->num_chan)   ? x->rows : x->num_chan;
  const int maxsamp = (x->cols < x->block_size) ? x->cols : x->block_size;
  t_atom *lptr=x->list_in;

  for (chan=0; chan<maxchan; chan++) {
    for (samp=0; samp<maxsamp; samp++) {
      x->sig_out[chan][samp]=atom_getfloat(&lptr[samp]);
    }
    lptr+=x->cols;
  }

  // zero missing signal samples
  lptr=x->list_in;
  for (chan=0; chan<maxchan; chan++) {
    for (; samp<x->block_size; samp++) {
      x->sig_out[chan][samp]=0;
      lptr+=x->cols;
    }
  }
  // zero missing channels
  for (chan=maxchan; chan<x->num_chan; chan++) {
    for (samp=0; samp<x->block_size; samp++) {
      x->sig_out[chan][samp]=0;
    }
    lptr+=x->cols;
  }

  // delete in the next dsp cycle, unless overwritten
  // by new matrix:
  x->perform_fcn=mtx_unpack_PerformSetInactive;

  return(arg+2);
}


void *mtx_unpack_new (t_floatarg f)
{
  int num_chan=1;
  mtx_unpack_tilde *x = (mtx_unpack_tilde*) pd_new(mtx_unpack_tilde_class);
  num_chan=(int)f;
  if ((num_chan<1) || (num_chan>MTX_PACK_MAXCHANNELS)) {
    num_chan=1;
  }
  x->num_chan=num_chan;
  x->sig_out=0;
  x->list_in=0;
  x->rows=0;
  x->cols=0;
  x->perform_fcn=mtx_unpack_PerformInactive;
  while (num_chan--) {
    outlet_new(&x->x_obj, &s_signal);
  }
  x->sig_out = (t_float**)getbytes(sizeof(t_float*)*x->num_chan);

  return (void *) x;
}
void mtx_unpack_delete (mtx_unpack_tilde *x)
{
  if (x->sig_out) {
    freebytes (x->sig_out, x->num_chan * sizeof (t_float));
  }
}
static void mtx_unpack_matrix (mtx_unpack_tilde *x, t_symbol *s,
                                  int argc, t_atom *argv)
{
  int rows, cols;
  x->rows=0;
  x->cols=0;
  if(iemmatrix_check(x, argc, argv, 0))return;
  rows=(int) atom_getfloat (argv++);
  cols=(int) atom_getfloat (argv++);
  argc-=2;
  x->rows=rows;
  x->cols=cols;
  x->list_in=argv;
  x->perform_fcn=mtx_unpack_PerformActive;
}

static void mtx_unpack_dsp (mtx_unpack_tilde *x, t_signal **sp)
{
  int chan;
  for (chan=0; chan<x->num_chan; chan++) {
    x->sig_out[chan]=sp[chan]->s_vec;
  }

  x->block_size=sp[0]->s_n;
  x->perform_fcn=mtx_unpack_PerformInactive;
  dsp_add(mtx_unpack_Perform,1,x);
}

void mtx_unpack_tilde_setup (void)
{
  mtx_unpack_tilde_class = class_new(gensym("mtx_unpack~"),
                                     (t_newmethod)mtx_unpack_new, (t_method) mtx_unpack_delete,
                                     sizeof(mtx_unpack_tilde), CLASS_DEFAULT, A_DEFFLOAT, 0);
  class_addmethod (mtx_unpack_tilde_class, (t_method) mtx_unpack_matrix,
                   gensym("matrix"),A_GIMME,0);
  class_addmethod (mtx_unpack_tilde_class, (t_method) mtx_unpack_dsp,
                   gensym("dsp"),0);
}

void iemtx_unpack__setup(void)
{
  mtx_unpack_tilde_setup();
}
