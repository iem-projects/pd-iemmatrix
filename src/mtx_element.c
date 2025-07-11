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


/* mtx_element */
static t_class *mtx_element_class;

static void mtx_element_list2(t_matrixobj *x, t_floatarg f1, t_floatarg f2)
{
  int r = f1, c= f2;
  if(r<0) {
    r=0;
  }
  if(c<0) {
    c=0;
  }
  x->current_row = r;
  x->current_col = c;
}
static void mtx_element_matrix(t_matrixobj *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  matrix_matrix2(x, &x->m, argc, argv);
  matrixobj_bang(x);
}
static void mtx_element_float(t_matrixobj *x, t_floatarg f)
{
  if(x->current_col>x->m.col || x->current_row>x->m.row) {
    pd_error(x,"[mtx_element]: element position exceeds matrix dimensions");
    return;
  }
  if(x->current_row == 0 && x->current_col == 0) {
    matrix_set(&x->m, f);
    matrixobj_bang(x);
    return;
  }
  if(x->current_row && x->current_col) {
    SETFLOAT(x->m.atombuffer+1+(x->current_row-1)*x->m.col+x->current_col, f);
  } else {
    t_atom *ap=x->m.atombuffer+2;
    int count;
    if (!x->current_col) {
      ap+=x->m.col*(x->current_row-1);
      count=x->m.col;
      while(count--) {
        SETFLOAT(&ap[count], f);
      }
    } else {
      ap+=x->current_col-1;
      count=x->m.row;
      while(count--) {
        SETFLOAT(&ap[count*x->m.col], f);
      }
    }
  }
  matrixobj_bang(x);
}

static void *mtx_element_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_element_class);
  int i, j, q;
  (void)s; /* unused */
  outlet_new(&x->x_obj, 0);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym(""));
  x->current_row=x->current_col=0;
  x->m.col=x->m.row=0;
  x->m.atombuffer=0;
  switch (argc) {
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
  case 4:
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
    q = atom_getfloat(argv++);
    if(q<0) {
      q=0;
    }
    x->current_row=q;
    q = atom_getfloat(argv++);
    if(q<0) {
      q=0;
    }
    x->current_col=q;
    break;
  default:
    ;
  }
  return (x);
}
void mtx_element_setup(void)
{
  mtx_element_class = class_new(gensym("mtx_element"),
                                (t_newmethod)mtx_element_new,
                                (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addbang  (mtx_element_class, matrixobj_bang);
  class_addfloat (mtx_element_class, mtx_element_float);
  class_addmethod(mtx_element_class, (t_method)mtx_element_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addmethod(mtx_element_class, (t_method)mtx_element_list2, gensym(""),
                  A_FLOAT, A_FLOAT, 0);

}


void iemtx_element_setup(void)
{
  mtx_element_setup();
}
