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

/* mtx_check */
static t_class *mtx_check_class;

static void mtx_check_matrix(t_matrixobj *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int row=atom_getfloat(argv);
  int col=atom_getfloat(argv+1);
  t_atom *ap;
  int length=row*col, n;
  (void)s; /* unused */

  argc-=2;

  if(length>argc) {
    /* sparse matrix */
    adjustsize(x, &x->m, row, col);
    matrix_set(&x->m, 0);
    argv+=2;
    ap=x->m.atombuffer+2;
    n=argc;
    while(n--) {
      t_float f = atom_getfloat(argv++);
      SETFLOAT(ap, f);
      ap++;
    }
    matrixobj_bang(x);
  } else {
    SETFLOAT(argv, row);
    SETFLOAT(argv+1, col);
    ap=argv+2;
    n=length;
    while(n--) {
      t_float f = atom_getfloat(ap);
      SETFLOAT(ap, f);
      ap++;
    }
    outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), length+2, argv);
  }
}

static void *mtx_check_new()
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_check_class);
  outlet_new(&x->x_obj, 0);
  return (x);
}

void mtx_check_setup(void)
{
  mtx_check_class = class_new(gensym("mtx_check"),
                              (t_newmethod)mtx_check_new,
                              (t_method)matrixobj_free, sizeof(t_matrixobj), 0, 0);
  class_addbang  (mtx_check_class, matrixobj_bang);
  class_addmethod(mtx_check_class, (t_method)mtx_check_matrix,
                  gensym("matrix"), A_GIMME, 0);

}


void iemtx_check_setup(void)
{
  mtx_check_setup();
}
