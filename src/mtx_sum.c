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

/* mtx_sum */
/* column-wise sum
 */
static t_class *mtx_sum_class;
static void mtx_sum_matrix(t_matrixobj *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  int row, col;
  int n;
  t_atom *ap = 0, *dummy=0;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);

  dummy = ap = (t_atom *)getbytes(col * sizeof(t_atom));

  for(n=0; n<col; n++, dummy++) {
    int i=row;
    t_float f=0.f;
    t_atom*ap2=argv+n;
    while(i--) {
      f+=atom_getfloat(ap2+col*i);
    }
    SETFLOAT(dummy, f);
  }

  outlet_list(x->x_obj.ob_outlet, gensym("prod"), col, ap);

  freebytes(ap, (col * sizeof(t_atom)));
}
static void mtx_sum_list(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  t_float f=0.f;
  (void)s; /* unused */
  while(argc--) {
    f+=atom_getfloat(argv++);
  }
  outlet_float(x->x_obj.ob_outlet, f);
}


static void *mtx_sum_new(void)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_sum_class);
  outlet_new(&x->x_obj, 0);
  x->m.row = x->m.col = 0;
  x->m.atombuffer   = 0;

  return (x);
}
void mtx_sum_setup(void)
{
  mtx_sum_class = class_new(gensym("mtx_sum"), (t_newmethod)mtx_sum_new,
                            (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addlist  (mtx_sum_class, mtx_sum_list);
  class_addmethod(mtx_sum_class, (t_method)mtx_sum_matrix, gensym("matrix"),
                  A_GIMME, 0);

}
void iemtx_sum_setup(void)
{
  mtx_sum_setup();
}
