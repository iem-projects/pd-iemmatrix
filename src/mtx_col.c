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

/* mtx_col */
static t_class *mtx_col_class;

static void mtx_col_float(t_matrixobj *x, t_floatarg f)
{
  int i = f;
  if(i<0) {
    i=0;
  }
  x->current_col = i;
}
static void mtx_col_matrix(t_matrixobj *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  matrix_matrix2(x, &x->m, argc, argv);
  matrixobj_bang(x);
}
static void mtx_col_list(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s; /* unused */
  if (argc==1) {
    t_float f=atom_getfloat(argv);
    t_atom *ap=x->m.atombuffer+1+x->current_col;
    if (x->current_col>x->m.col) {
      pd_error(x, "[mtx_col]: too high a column is to be set");
      return;
    }
    if (x->current_col) {
      int n=x->m.row;
      while(n--) {
        SETFLOAT(ap, f);
        ap+=x->m.row+1;
      }
    }
    matrixobj_bang(x);
    return;
  }

  if (argc<x->m.row) {
    pd_error(x, "[mtx_col]: column length is too small for %dx%d-matrix", x->m.row,
         x->m.col);
    return;
  }
  if (x->current_col>x->m.col) {
    pd_error(x, "[mtx_col]: too high a column is to be set");
    return;
  }
  if(x->current_col) {
    int r=x->m.row;
    t_atom *ap=x->m.atombuffer+1+x->current_col;
    while(r--) {
      SETFLOAT(&ap[(x->m.row-r-1)*x->m.col], atom_getfloat(argv++));
    }
  }  else {
    int r=x->m.row;
    t_atom *ap=x->m.atombuffer+2;
    while (r--) {
      t_float f=atom_getfloat(argv++);
      int c=x->m.col;
      while(c--) {
        SETFLOAT(ap, f);
        ap++;
      }
    }
  }
  matrixobj_bang(x);
}
static void *mtx_col_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_col_class);
  int i, j, q;
  (void)s; /* unused */
  outlet_new(&x->x_obj, 0);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym(""));
  x->current_col=0;
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
    x->current_col=q;
  }
  return (x);
}
void mtx_col_setup(void)
{
  mtx_col_class = class_new(gensym("mtx_col"), (t_newmethod)mtx_col_new,
                            (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addbang  (mtx_col_class, matrixobj_bang);
  class_addlist  (mtx_col_class, mtx_col_list);
  class_addmethod(mtx_col_class, (t_method)mtx_col_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addmethod(mtx_col_class, (t_method)mtx_col_float, gensym(""),
                  A_FLOAT, 0);

}

void iemtx_col_setup(void)
{
  mtx_col_setup();
}
