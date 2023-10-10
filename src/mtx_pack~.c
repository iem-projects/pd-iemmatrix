#include "iemmatrix.h"

#ifndef CLASS_MULTICHANNEL
# define CLASS_MULTICHANNEL 0
#endif
#define MTX_PACK_MAXCHANNELS 200

static t_class *mtx_pack_tilde_class;

typedef struct _mtx_pack_tilde {
  t_object x_obj;
  int block_size;
  size_t num_ports; /* number of inlet~s (0 for 1 multichannel signal) */
  size_t num_channels; /* number of input signals */
  t_sample **sig_in;
  t_atom *list_out;
  t_outlet *message_outlet;

  t_clock*clock;

} mtx_pack_tilde;

static void mTxPackTildeOut(mtx_pack_tilde*x) {
  outlet_anything(x->message_outlet,gensym("matrix"),
                  x->block_size*x->num_channels+2,x->list_out);
}
void *newMtxPackTilde (t_floatarg f)
{
  int deferred = 0;
  int num_ports=1;
  mtx_pack_tilde *x = (mtx_pack_tilde*) pd_new(mtx_pack_tilde_class);
  num_ports=(int)f;
  if ((num_ports < 1) || (num_ports>MTX_PACK_MAXCHANNELS)) {
    num_ports=1;
  }
  x->num_ports=num_ports;
  while (num_ports--) {
    signalinlet_new(&x->x_obj, 0);
  }
  x->message_outlet=(t_outlet*)outlet_new(&x->x_obj,&s_list);

  if(deferred)
    x->clock = clock_new(x, (t_method)mTxPackTildeOut);

  return (void *) x;
}
void deleteMtxPackTilde (mtx_pack_tilde *x)
{
  if (x->sig_in) {
    freebytes (x->sig_in, x->num_channels * sizeof (*x->sig_in));
  }
  if (x->list_out) {
    freebytes (x->list_out, (x->num_channels * x->block_size + 2)*sizeof(*x->list_out));
  }
  if(x->clock)
    clock_free(x->clock);
}
static t_int *mTxPackTildePerform (t_int *arg)
{
  mtx_pack_tilde *x = (mtx_pack_tilde *) (arg[1]);
  int chan;
  int samp;
  t_atom *lptr=x->list_out+2;
  SETFLOAT(x->list_out,(t_float)x->num_channels);
  SETFLOAT(x->list_out+1,(t_float)x->block_size);

  for (chan=0; chan<x->num_channels; chan++) {
    for (samp=0; samp<x->block_size; samp++,lptr++) {
      SETFLOAT(lptr, x->sig_in[chan][samp]);
    }
  }

  if(x->clock)
    clock_delay(x->clock, 0);
  else
    mTxPackTildeOut(x);

  return(arg+2);

}

static void mTxPackTildeDsp (mtx_pack_tilde *x, t_signal **sp)
{
  size_t i, chan = x->num_ports;
  int block_size=sp[0]->s_n;

  if(x->sig_in)
    freebytes(x->sig_in, sizeof(*x->sig_in) * x->num_channels);
  x->sig_in = 0;

#if CLASS_MULTICHANNEL
  /* with multichannels, we add concatenate all channels from all ports */
  chan = 0;
  for(i=0; i<x->num_ports; i++) {
    chan += sp[i]->s_nchans;
  }
#endif

  x->num_channels = chan;
  x->sig_in = (t_sample**)getbytes(sizeof(*x->sig_in)*x->num_channels);

#if CLASS_MULTICHANNEL
  for(i=0, chan=0; i<x->num_ports; i++) {
    int j;
    for(j=0; j<sp[i]->s_nchans; j++) {
      x->sig_in[chan++] = sp[i]->s_vec + j * block_size;
    }
  }
#else
  for (i=0; i < chan; i++) {
    x->sig_in[i]=sp[i]->s_vec;
  }
#endif

  x->block_size=block_size;
  x->list_out = (t_atom*) getbytes ((x->num_channels * x->block_size + 2)
                                    * sizeof(*x->list_out));
  dsp_add(mTxPackTildePerform,1,x);
}

void mtx_pack_tilde_setup (void)
{
  int flags = CLASS_NOINLET;
  if (CLASS_MULTICHANNEL && iemmatrix_getpdfun("signal_setmultiout")) {
    flags |= CLASS_MULTICHANNEL;
  }

  mtx_pack_tilde_class = class_new(gensym("mtx_pack~"),
                                   (t_newmethod)newMtxPackTilde, (t_method)deleteMtxPackTilde,
                                   sizeof(mtx_pack_tilde),
				   flags,
				   A_DEFFLOAT, 0);
  class_addmethod (mtx_pack_tilde_class, (t_method) mTxPackTildeDsp,
                   gensym("dsp"),0);
}

void iemtx_pack__setup(void)
{
  mtx_pack_tilde_setup();
}
