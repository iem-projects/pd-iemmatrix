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

/* mtx_rmstodb: B=log(A); B[n,m]=e^A[n,m]  */

static t_class *mtx_rmstodb_class;

static void mtx_rmstodb_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  int row, col;
  int n = argc-2;
  t_atom *m;
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);

  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=(f<0)?0.:(100+20./LOGTEN * log(f));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_rmstodb_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=(f<0)?0.:(100+20./LOGTEN * log(f));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_rmstodb_new()
{
  /* element log */
  t_matrix *x = (t_matrix *)pd_new(mtx_rmstodb_class);
  outlet_new(&x->x_obj, 0);
  return(x);
}

void mtx_rmstodb_setup(void)
{
  mtx_rmstodb_class = class_new(gensym("mtx_rmstodb"),
                                (t_newmethod)mtx_rmstodb_new, (t_method)mtx_binmtx_free,
                                sizeof(t_mtx_binmtx), 0, 0);
  class_addmethod(mtx_rmstodb_class, (t_method)mtx_rmstodb_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_rmstodb_class, mtx_rmstodb_list);
  class_addbang  (mtx_rmstodb_class, mtx_binmtx_bang);


}

void iemtx_rmstodb_setup(void)
{
  mtx_rmstodb_setup();
}
