/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */
#include "iemmatrix.h"


/* mtx_mean */
static t_class *mtx_mean_class;

static void mtx_mean_matrix(t_matrixobj *x, t_symbol *s, int argc,
                            t_atom *argv)
{
  if(argc<2)return;
  int row=atom_getfloat(argv++);
  int col=atom_getfloat(argv++);
  t_atom *ip, *op;
  int c=col, r;
  t_float sum;
  t_float factor=1./row;
  (void)s; /* unused */
  adjustsize(x, &x->m, 1, col);
  op=x->m.atombuffer;

  while(c--) {
    sum=0;
    ip=argv+col-c-1;
    r=row;
    while(r--) {
      sum+=atom_getfloat(ip+col*r);
    }
    SETFLOAT(op, sum*factor);
    op++;
  }
  outlet_list(x->x_obj.ob_outlet, gensym("row"), col, x->m.atombuffer);
}

static void *mtx_mean_new(void)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_mean_class);
  outlet_new(&x->x_obj, 0);
  return (x);
}
void mtx_mean_setup(void)
{
  mtx_mean_class = class_new(gensym("mtx_mean"), (t_newmethod)mtx_mean_new,
                             (t_method)matrixobj_free, sizeof(t_matrixobj), 0, 0, 0);
  class_addmethod(mtx_mean_class, (t_method)mtx_mean_matrix,
                  gensym("matrix"), A_GIMME, 0);

}

void iemtx_mean_setup(void)
{
  mtx_mean_setup();
}
