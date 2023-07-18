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

/* mtx_not: B=!A; */

#define MTX_ALMOSTZERO 1e-19

static t_class *mtx_not_class;

static void mtx_not_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  int row, col, n;
  t_atom *m;
  if(iemmatrix_check(x, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  n = row*col;

  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f = atom_getfloat(argv++);
    SETFLOAT(m, (t_float)(f<MTX_ALMOSTZERO&&f>-MTX_ALMOSTZERO));
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_not_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                         t_atom *argv)
{
  int n=argc;
  t_atom *m;

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    t_float f = atom_getfloat(argv++);
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = (t_float)(f<MTX_ALMOSTZERO&&f>-MTX_ALMOSTZERO);
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_not_new(t_symbol *s)
{
  /* element not */
  t_matrix *x = (t_matrix *)pd_new(mtx_not_class);
  outlet_new(&x->x_obj, 0);
  x->col = x->row = 0;
  x->atombuffer = 0;
  return(x);
}

void mtx_not_setup(void)
{
  mtx_not_class = class_new(gensym("mtx_not"), (t_newmethod)mtx_not_new,
                            (t_method)mtx_binmtx_free,
                            sizeof(t_mtx_binmtx), 0, A_GIMME, 0);
  class_addcreator((t_newmethod)mtx_not_new, gensym("mtx_!"), A_GIMME,0);
  class_addmethod(mtx_not_class, (t_method)mtx_not_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addlist  (mtx_not_class, mtx_not_list);
  class_addbang  (mtx_not_class, mtx_binmtx_bang);


}

void iemtx_not_setup(void)
{
  mtx_not_setup();
}
