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

/* mtx_abs: B=abs(A); */

static t_class *mtx_abs_class;

static void mtx_abs_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  t_atom *m;
  int n;
  int row, col;
  (void)s; /* unused */

  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  n = row*col;
  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f = atom_getfloat(argv++);
    SETFLOAT(m, (t_float)fabs(f));
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_abs_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                         t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = (t_float)fabs(atom_getfloat(argv++));
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_abs_new(t_symbol *s)
{
  /* element abs */
  t_matrix *x = (t_matrix *)pd_new(mtx_abs_class);
  (void)s; /* unused */

  outlet_new(&x->x_obj, 0);
  return(x);
}

void mtx_abs_setup(void)
{
  mtx_abs_class = class_new(gensym("mtx_abs"), (t_newmethod)mtx_abs_new,
                            (t_method)mtx_binmtx_free,
                            sizeof(t_mtx_binmtx), 0, A_GIMME, 0);
  class_addmethod(mtx_abs_class, (t_method)mtx_abs_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addlist  (mtx_abs_class, mtx_abs_list);
  class_addbang  (mtx_abs_class, mtx_binmtx_bang);


}

void iemtx_abs_setup(void)
{
  mtx_abs_setup();
}
