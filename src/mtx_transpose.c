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

/* mtx_transpose */
static t_class *mtx_transpose_class;

static void mtx_transpose_matrix(t_matrixobj *x, t_symbol *s, int argc,
                                 t_atom *argv)
{
  int row, col;
  t_atom *ap;
  int r, c;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv+0);
  col=atom_getint(argv+1);
  if (col*row!=x->m.col*x->m.row) {
    freebytes(x->m.atombuffer, (x->m.col*x->m.row+2)*sizeof(t_atom));
    x->m.atombuffer = (t_atom *)getbytes((row*col+2)*sizeof(t_atom));
  }
  ap = x->m.atombuffer+2;
  setdimen(&x->m, col, row);
  r = row;
  while(r--) {
    c=col;
    while(c--) {
      t_float f = atom_getfloat(argv+2+r*col+c);
      SETFLOAT(ap+c*row+r, f);
    }
  }

  matrixobj_bang(x);
}

static void *mtx_transpose_new()
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_transpose_class);
  outlet_new(&x->x_obj, 0);
  return (x);
}
void mtx_transpose_setup(void)
{
  mtx_transpose_class = class_new(gensym("mtx_transpose"),
                                  (t_newmethod)mtx_transpose_new,
                                  (t_method)matrixobj_free, sizeof(t_matrixobj), 0, 0);
  class_addbang  (mtx_transpose_class, matrixobj_bang);
  class_addmethod(mtx_transpose_class, (t_method)mtx_transpose_matrix,
                  gensym("matrix"), A_GIMME, 0);

}

void iemtx_transpose_setup(void)
{
  mtx_transpose_setup();
}
