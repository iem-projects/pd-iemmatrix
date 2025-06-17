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

/* mtx_int: B=int(A); */

static t_class *mtx_int_class;

static void mtx_int_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  int row, col, n;
  t_atom *m;
  (void)s; /* unused */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  n=row*col;

  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f = (t_float)atom_getint(argv++);
    SETFLOAT(m, f);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_int_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                         t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = (t_float)atom_getint(argv++);
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_int_new()
{
  /* element int */
  t_matrix *x = (t_matrix *)pd_new(mtx_int_class);
  outlet_new(&x->x_obj, 0);
  return(x);
}

void mtx_int_setup(void)
{
  mtx_int_class = class_new(gensym("mtx_int"), (t_newmethod)mtx_int_new,
                            (t_method)mtx_binmtx_free,
                            sizeof(t_mtx_binmtx), 0, 0);
  class_addmethod(mtx_int_class, (t_method)mtx_int_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addlist  (mtx_int_class, mtx_int_list);
  class_addbang  (mtx_int_class, mtx_binmtx_bang);


}

void iemtx_int_setup(void)
{
  mtx_int_setup();
}
