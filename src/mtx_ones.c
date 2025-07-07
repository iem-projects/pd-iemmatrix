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

/* mtx_ones */
static t_class *mtx_ones_class;
static void *mtx_ones_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_ones_class);
  int col=0, row=0;
  (void)s; /* unused */
  outlet_new(&x->x_obj, 0);
  x->m.row = x->m.col = 0;
  x->m.atombuffer   = 0;
  switch(argc) {
  case 0:
    break;
  case 1:
    col=row=atom_getfloat(argv);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
  }
  if(col<0) {
    col=0;
  }
  if(row<0) {
    row=0;
  }
  if (col && row) {
    x->m.atombuffer = (t_atom *)getbytes((col*row+2)*sizeof(t_atom));
    setdimen(&x->m, row, col);
    matrix_set(&x->m, 1);
  }
  return (x);
}
void mtx_ones_setup(void)
{
  mtx_ones_class = class_new(gensym("mtx_ones"), (t_newmethod)mtx_ones_new,
                             (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addlist(mtx_ones_class, matrixobj_ones);
  class_addbang(mtx_ones_class, matrixobj_bang);
  class_addmethod(mtx_ones_class, (t_method)matrixobj_ones, gensym("matrix"),
                  A_GIMME, 0);


}
void iemtx_ones_setup(void)
{
  mtx_ones_setup();
}
