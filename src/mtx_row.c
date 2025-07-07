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

/* mtx_row */
static t_class *mtx_row_class;

static void mtx_row_float(t_matrixobj *x, t_floatarg f)
{
  int i = f;
  if(i<0) {
    i=0;
  }
  x->current_row = i;
}
static void mtx_row_matrix(t_matrixobj *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  matrix_matrix2(x, &x->m, argc, argv);
  matrixobj_bang(x);
}
static void mtx_row_list(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s; /* unused */
  if (argc==1) {
    t_float f=atom_getfloat(argv);
    t_atom *ap=x->m.atombuffer+2+(x->current_row-1)*x->m.col;
    if (x->current_row>x->m.row) {
      pd_error(x, "[mtx_row]: too high a row is to be set");
      return;
    }
    if (x->current_row) {
      int n=x->m.col;
      while(n--) {
        SETFLOAT(ap, f);
        ap++;
      }
    }
    matrixobj_bang(x);
    return;
  }

  if (argc<x->m.col) {
    pd_error(x, "[mtx_row]: row length is too small for %dx%d-matrix", x->m.row, x->m.col);
    return;
  }
  if (x->current_row>x->m.row) {
    pd_error(x, "[mtx_row]: too high a row is to be set");
    return;
  }
  if(x->current_row) {
    memcpy(x->m.atombuffer+2+(x->current_row-1)*x->m.col, argv,
           x->m.col*sizeof(t_atom));
  }  else {
    int r=x->m.row;
    while(r--) {
      memcpy(x->m.atombuffer+2+r*x->m.col, argv, x->m.col*sizeof(t_atom));
    }
  }
  matrixobj_bang(x);
}
static void *mtx_row_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_row_class);
  int i, j, q;
  (void)s; /* unused */

  outlet_new(&x->x_obj, 0);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym(""));
  x->current_row=0;
  x->m.col=x->m.row=0;
  x->m.atombuffer=0;
  switch (argc) {
  case 0:
    break;
  case 1:
    i = atom_getfloat(argv);
    if (i<0) {
      i=0;
    }
    if(i) {
      adjustsize(x, &x->m, i, i);
    }
    matrix_set(&x->m, 0);
    break;
  case 2:
    i = atom_getfloat(argv++);
    if(i<0) {
      i=0;
    }
    j = atom_getfloat(argv++);
    if(j<0) {
      j=0;
    }
    if(i && j) {
      adjustsize(x, &x->m, i, j);
    }
    matrix_set(&x->m, 0);
    break;
  default:
    i = atom_getfloat(argv++);
    if(i<0) {
      i=0;
    }
    j = atom_getfloat(argv++);
    if(j<0) {
      j=0;
    }
    q = atom_getfloat(argv++);
    if(q<0) {
      q=0;
    }
    if(i && j) {
      adjustsize(x, &x->m, i, j);
    }
    matrix_set(&x->m, 0);
    x->current_row=q;
  }
  return (x);
}
void mtx_row_setup(void)
{
  mtx_row_class = class_new(gensym("mtx_row"), (t_newmethod)mtx_row_new,
                            (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addbang  (mtx_row_class, matrixobj_bang);
  class_addlist  (mtx_row_class, mtx_row_list);
  class_addmethod(mtx_row_class, (t_method)mtx_row_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addmethod(mtx_row_class, (t_method)mtx_row_float, gensym(""),
                  A_FLOAT, 0);

}


void iemtx_row_setup(void)
{
  mtx_row_setup();
}
