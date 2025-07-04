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

#define LOGTEN 2.302585092994

/* mtx_powtodb: B=log(A); B[n,m]=e^A[n,m]  */

static t_class *mtx_powtodb_class;

static void mtx_powtodb_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  int row, col, n;
  t_atom *m;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  n = row*col;

  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=(f<0)?0.:(100+10./LOGTEN * log(f));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_powtodb_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=(f<0)?0.:(100+10./LOGTEN * log(f));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_powtodb_new()
{
  /* element log */
  t_matrix *x = (t_matrix *)pd_new(mtx_powtodb_class);
  outlet_new(&x->x_obj, 0);
  return(x);
}

void mtx_powtodb_setup(void)
{
  mtx_powtodb_class = class_new(gensym("mtx_powtodb"),
                                (t_newmethod)mtx_powtodb_new, (t_method)mtx_binmtx_free,
                                sizeof(t_mtx_binmtx), 0, 0);
  class_addmethod(mtx_powtodb_class, (t_method)mtx_powtodb_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_powtodb_class, mtx_powtodb_list);
  class_addbang  (mtx_powtodb_class, mtx_binmtx_bang);


}

void iemtx_powtodb_setup(void)
{
  mtx_powtodb_setup();
}
