/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly refering to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */
#include "iemmatrix.h"

/* mtx_cos: B=cos(A); */

static t_class *mtx_cos_class;

static void mtx_cos_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  int row, col;
  t_atom *m;
  int n = argc-2;
  if(iemmatrix_check(x, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f = (t_float)cos(atom_getfloat(argv++));
    SETFLOAT(m, f);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_cos_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                         t_atom *argv)
{
  int n=argc;
  t_atom *m;

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = (t_float)cos(atom_getfloat(argv++));
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_cos_new(t_symbol *s)
{
  /* element cos */
  t_matrix *x = (t_matrix *)pd_new(mtx_cos_class);
  outlet_new(&x->x_obj, 0);
  x->col = x->row = 0;
  x->atombuffer = 0;
  return(x);
}

void mtx_cos_setup(void)
{
  mtx_cos_class = class_new(gensym("mtx_cos"), (t_newmethod)mtx_cos_new,
                            (t_method)mtx_binmtx_free,
                            sizeof(t_mtx_binmtx), 0, A_GIMME, 0);
  class_addmethod(mtx_cos_class, (t_method)mtx_cos_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addlist  (mtx_cos_class, mtx_cos_list);
  class_addbang  (mtx_cos_class, mtx_binmtx_bang);


}

void iemtx_cos_setup(void)
{
  mtx_cos_setup();
}
