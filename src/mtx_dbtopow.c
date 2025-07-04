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

/* mtx_dbtopow: B=log(A); B[n,m]=e^A[n,m]  */

static t_class *mtx_dbtopow_class;

static void mtx_dbtopow_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  int row, col;
  t_atom *m;
  int n = argc-2;
  (void)s; /* unused */
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);

  adjustsize(&x->m, row, col);
  m =  x->m.atombuffer+2;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=0;
    f=(f>485)?485:f;
    v=(f<=0)?0:exp((LOGTEN*0.1) * (f-100.));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}

static void mtx_dbtopow_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */

  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    t_float f=atom_getfloat(argv++);
    t_float v=0;
    f=(f>485)?485:f;
    v=(f<=0)?0:exp((LOGTEN*0.1) * (f-100.));
    SETFLOAT(m, (v<0)?0:v);
    m++;
  }

  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_dbtopow_new(void)
{
  /* element log */
  t_matrix *x = (t_matrix *)pd_new(mtx_dbtopow_class);
  outlet_new(&x->x_obj, 0);
  return(x);
}

void mtx_dbtopow_setup(void)
{
  mtx_dbtopow_class = class_new(gensym("mtx_dbtopow"),
                                (t_newmethod)mtx_dbtopow_new, (t_method)mtx_binmtx_free,
                                sizeof(t_mtx_binmtx), 0, 0);
  class_addmethod(mtx_dbtopow_class, (t_method)mtx_dbtopow_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_dbtopow_class, mtx_dbtopow_list);
  class_addbang  (mtx_dbtopow_class, mtx_binmtx_bang);


}

void iemtx_dbtopow_setup(void)
{
  mtx_dbtopow_setup();
}
