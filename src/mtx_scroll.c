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

/* mtx_scroll */
/* scroll the rows */
static t_class *mtx_scroll_class;

static void mtx_scroll_matrix(t_mtx_binscalar *x, t_symbol *s, int argc,
                              t_atom *argv)
{
  int row, col, rowscroll;
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  row=atom_getfloat(argv++);
  col=atom_getfloat(argv++);
  rowscroll = ((int)x->f%row+row)%row;
  adjustsize(x, &x->m, row, col);

  memcpy(x->m.atombuffer+2, argv+(row-rowscroll)*col,
         rowscroll*col*sizeof(t_atom));
  memcpy(x->m.atombuffer+2+rowscroll*col, argv,
         (row-rowscroll)*col*sizeof(t_atom));

  matrixobj_bang((t_matrixobj*)x);
}

static void *mtx_scroll_new(t_floatarg f)
{
  t_mtx_binscalar *x = (t_mtx_binscalar *)pd_new(mtx_scroll_class);
  floatinlet_new(&x->x_obj, &(x->f));
  outlet_new(&x->x_obj, 0);

  x->f=f;
  return (x);
}
void mtx_scroll_setup(void)
{
  mtx_scroll_class = class_new(gensym("mtx_scroll"),
                               (t_newmethod)mtx_scroll_new,
                               (t_method)matrixobj_free, sizeof(t_mtx_binscalar), 0, A_DEFFLOAT, 0);
  class_addbang  (mtx_scroll_class, matrixobj_bang);
  class_addmethod(mtx_scroll_class, (t_method)mtx_scroll_matrix,
                  gensym("matrix"), A_GIMME, 0);

}
void iemtx_scroll_setup(void)
{
  mtx_scroll_setup();
}
