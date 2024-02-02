/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) 2008, Franz Zotter, IEM KUG Graz Austria
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute, IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"
#define MTX_PACK_MAXCHANNELS 200

#ifndef CLASS_MULTICHANNEL
# define CLASS_MULTICHANNEL 0
#endif
typedef void (*setmultiout_f)(t_signal **sig, int nchans);

static t_class *mtx_unpack_tilde_class;
static t_class *mtx_unpack_tilde_proxy;

typedef struct _proxy {
  t_object p_obj;
  struct _mtx_unpack_tilde*p_owner;
} t_proxy;



typedef struct _mtx_unpack_tilde {
  t_object x_obj;
  t_proxy *x_proxy;
  int x_dsp; /* whether DSP is currently running */
  int rows, cols;
  int block_size;
  unsigned int num_chan, num_ports;
  t_float **sig_out;
  t_atom *list_in;
  t_int *(*perform_fcn)(t_int*);

  /* when doing multichannel, this is Pd>=0.54's signal_setmultiout(); otherwise NULL */
  setmultiout_f x_setmultiout;
} mtx_unpack_tilde;

static void proxy_dspstopped(t_proxy*p) {
  p->p_owner->x_dsp = 0;
}

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


void *mtx_unpack_new (t_symbol*s, int argc, t_atom*argv)
{
  setmultiout_f setmultiout = (CLASS_MULTICHANNEL)?iemmatrix_getpdfun("signal_setmultiout"):0;
  int want_multi = 0;
  int num_chan=(int)1;
  mtx_unpack_tilde *x = 0;  

  /* args:
     - (none): single-signal, n=1
     - <int> : single-signal, n=N
     - '-m': multi-signal, n=?
     - '-m <int>: multi-signal, n=N
  */
  if(argc && A_SYMBOL == argv->a_type) {
    if (atom_getsymbol(argv) == gensym("-m")) {
      want_multi = 1;
      argc--;
      argv++;
    } else {
      goto usage;
    }
  }
  if(argc && A_FLOAT != argv->a_type) {
    goto usage;
  }

  if(argv)
    num_chan = (int)atom_getfloat(argv);

  x = (mtx_unpack_tilde*) pd_new(mtx_unpack_tilde_class);

  /* DSP stopped proxy */
  x->x_proxy = (t_proxy*)pd_new(mtx_unpack_tilde_proxy);
  x->x_proxy->p_owner = x;
  pd_bind(&x->x_proxy->p_obj.ob_pd, gensym("pd-dsp-stopped"));

  /* check multichannel support */
  if(want_multi) {
    static int warn_multichannel = 1;
    x->x_setmultiout = setmultiout;
    if (warn_multichannel && !setmultiout) {
      if(CLASS_MULTICHANNEL) {
        int major, minor, bugfix;
        sys_getversion(&major, &minor, &bugfix);
        pd_error(x, "[%s] multichannel requested, but iemmatrix is running in Pd-%d.%d-%d, which doesn't support it", s->s_name, major, minor, bugfix);
      } else {
        pd_error(x, "[%s] multichannel requested, but iemmatrix was built against Pd-%d.%d-%d, which doesn't support it", s->s_name, PD_MAJOR_VERSION, PD_MINOR_VERSION, PD_BUGFIX_VERSION);
      }
    }
    warn_multichannel = 0;
  }
  if ((num_chan<1) || (num_chan>MTX_PACK_MAXCHANNELS)) {
    if(!want_multi)
      pd_error(x, "[mtx_unpack~] invalid number of channels (%d), default to 1.", num_chan);
    num_chan=1;
  }

  x->perform_fcn=mtx_unpack_PerformInactive;

  x->rows = num_chan;
  x->num_chan= ((want_multi && !x->x_setmultiout)?1:num_chan);
  x->num_ports = (want_multi?1:num_chan);
  
  x->sig_out = (t_float**)getbytes(sizeof(t_float*)*x->num_chan);

  num_chan = x->num_ports;
  while (num_chan--) {
    outlet_new(&x->x_obj, &s_signal);
  }

  return (void *) x;

 usage:
  pd_error(0, "[mtx_unpack~] bad arguments, use '<int:channels>' or '-m <int:channel>'");
  return 0;
}
void mtx_unpack_delete (mtx_unpack_tilde *x)
{
  if (x->sig_out) {
    freebytes (x->sig_out, x->num_chan * sizeof (t_float));
  }
  if(x->x_proxy) {
    pd_unbind(&x->x_proxy->p_obj.ob_pd, gensym("pd-dsp-stopped"));
    pd_free(&x->x_proxy->p_obj.ob_pd);
  }
}
static void mtx_unpack_matrix (mtx_unpack_tilde *x, t_symbol *s,
                                  int argc, t_atom *argv)
{
  int rows, cols;
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
  x->x_dsp = 1;
  x->block_size=sp[0]->s_n;

  if(x->x_setmultiout) {
#if CLASS_MULTICHANNEL
    if(x->rows != x->num_chan) {
      if (x->sig_out) {
        freebytes (x->sig_out, x->num_chan * sizeof (t_float));
      }
      x->num_chan = x->rows;
      if(x->rows < 1)
        x->num_chan = 1;
      x->sig_out = (t_float**)getbytes(sizeof(t_float*)*x->num_chan);
    }
    /* create multichannel output */
    x->x_setmultiout(&sp[0], x->num_chan);
    for(chan=0; chan<x->num_chan; chan++) {
	x->sig_out[chan] = sp[0]->s_vec + chan * x->block_size;
      }
#else
    pd_error(x, "BUG: multichannel enabled but not compile-time multichannel support");
    return;
#endif
  } else {
    for (chan=0; chan<x->num_chan; chan++) {
      x->sig_out[chan]=sp[chan]->s_vec;
    }
  }

  x->perform_fcn=mtx_unpack_PerformInactive;
  dsp_add(mtx_unpack_Perform, 1, x);
}


void mtx_unpack_tilde_setup (void)
{
  mtx_unpack_tilde_class = class_new(gensym("mtx_unpack~"),
                                     (t_newmethod)mtx_unpack_new, (t_method) mtx_unpack_delete,
                                     sizeof(mtx_unpack_tilde), CLASS_DEFAULT, A_GIMME, 0);
  class_addmethod (mtx_unpack_tilde_class, (t_method) mtx_unpack_matrix,
                   gensym("matrix"),A_GIMME,0);
  class_addmethod (mtx_unpack_tilde_class, (t_method) mtx_unpack_dsp,
                   gensym("dsp"),0);

  mtx_unpack_tilde_proxy = class_new(gensym("mtx_unpack~ proxy"),
                                     0, 0,
                                     sizeof(t_proxy),
                                     CLASS_PD,
                                     0);
  class_addbang(mtx_unpack_tilde_proxy, (t_method)proxy_dspstopped);
}

void iemtx_unpack__setup(void)
{
  mtx_unpack_tilde_setup();
}
