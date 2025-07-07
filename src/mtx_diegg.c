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
static t_class *mtx_diegg_class;
static void mtx_diegg_matrix(t_matrixobj *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int row, col, length, n;
  t_atom *ap = 0, *dummy=0;
  (void)s; /* unused */

  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getfloat(argv++);
  col=atom_getfloat(argv++);
  length=(col<row)?col:row;
  n=length;

  ap=(t_atom *)getbytes(length * sizeof(t_atom));
  dummy=ap;

  for(n=0; n<length; n++, dummy++) {
    int index=(n+1)*(col-1);
    SETFLOAT(dummy, atom_getfloat(argv+index));
  }
  outlet_list(x->x_obj.ob_outlet, gensym("diegg"), length, ap);

  freebytes(ap, (length * sizeof(t_atom)));
}
static void *mtx_diegg_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_diegg_class);
  outlet_new(&x->x_obj, 0);
  x->m.row = x->m.col = 0;
  x->m.atombuffer   = 0;

  if(!argc) {
    (void)s;
    return(x);
  }

  matrix_diegg(&x->x_obj, &x->m, argc, argv);

  return (x);
}
void mtx_diegg_setup(void)
{
  mtx_diegg_class = class_new(gensym("mtx_diegg"),
                              (t_newmethod)mtx_diegg_new,
                              (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addlist  (mtx_diegg_class, matrixobj_diegg);
  class_addbang  (mtx_diegg_class, matrixobj_bang);
  class_addmethod(mtx_diegg_class, (t_method)mtx_diegg_matrix,
                  gensym("matrix"), A_GIMME, 0);

}
void iemtx_diegg_setup(void)
{
  mtx_diegg_setup();
}
