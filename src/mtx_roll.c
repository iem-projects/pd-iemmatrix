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

/* mtx_roll */
/* roll the rows */
static t_class *mtx_roll_class;

static void mtx_roll_matrix(t_mtx_binscalar *x, t_symbol *s, int argc,
                            t_atom *argv)
{
  int row, col, colroll;
  t_atom *ap;
  int c;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  colroll = ((int)x->f%col+col)%col;

  adjustsize(x, &x->m, row, col);
  ap = x->m.atombuffer+2;

  c=col;
  while(c--) {
    t_atom *in  = argv+col-c-1;
    t_atom *out = ap  +(col-c-1+colroll)%col;
    int r = row;
    while (r--) {
      SETFLOAT(out, atom_getfloat(in));
      out+=col;
      in+=col;
    }

  }

  matrixobj_bang((t_matrixobj*)x);
}

static void *mtx_roll_new(t_floatarg f)
{
  t_mtx_binscalar *x = (t_mtx_binscalar *)pd_new(mtx_roll_class);
  floatinlet_new(&x->x_obj, &(x->f));
  outlet_new(&x->x_obj, 0);

  x->f=f;
  return (x);
}
void mtx_roll_setup(void)
{
  mtx_roll_class = class_new(gensym("mtx_roll"), (t_newmethod)mtx_roll_new,
                             (t_method)matrixobj_free, sizeof(t_mtx_binscalar), 0, A_DEFFLOAT, 0);
  class_addbang  (mtx_roll_class, matrixobj_bang);
  class_addmethod(mtx_roll_class, (t_method)mtx_roll_matrix,
                  gensym("matrix"), A_GIMME, 0);

}
void iemtx_roll_setup(void)
{
  mtx_roll_setup();
}
