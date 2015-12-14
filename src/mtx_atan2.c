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

/*
  mtx_atan2
*/

/* mtx_atan2 */
static t_class *mtx_atan2_class, *mtx_atan2scalar_class;

static void mtx_atan2scalar_matrix(t_mtx_binscalar *x, t_symbol *s,
                                   int argc, t_atom *argv)
{
  int n=argc-2;
  int row=atom_getfloat(argv), col=atom_getfloat(argv+1);

  t_float offset=x->f;
  t_atom *buf;
  t_atom *ap=argv+2;
  if(iemmatrix_check(x, argc, argv, IEMMATRIX_CHECK_CRIPPLED))return;
  adjustsize(&x->m, row, col);

  buf=x->m.atombuffer+2;

  while(n--) {
    buf->a_type = A_FLOAT;
    buf++->a_w.w_float = atan2f(atom_getfloat(ap),offset);
    ap++;
  }
  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}
static void mtx_atan2scalar_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                                 t_atom *argv)
{
  int n=argc;
  t_atom *m;
  t_float offset = x->f;
  adjustsize(&x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = atan2f(atom_getfloat(argv),offset);
    argv++;
  }
  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void mtx_atan2_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int row=atom_getfloat(argv);
  int col=atom_getfloat(argv+1);
  t_atom *m;
  t_atom *m1 = argv+2;
  t_atom *m2 = x->m2.atombuffer+2;
  int n = argc-2;
  if(iemmatrix_check(x, argc, argc, 0))return;

  if (!(x->m2.col*x->m2.row)) {
    outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc, argv);
    return;
  }

  if ((col!=x->m2.col)||(row!=x->m2.row)) {
    post("mtx_atan2: matrix dimensions do not match");
    /* LATER SOLVE THIS */
    return;
  }
  adjustsize(&x->m, row, col);
  m = x->m.atombuffer+2;

  while(n--) {
    t_float f1=atom_getfloat(m1++);
    t_float f2=atom_getfloat(m2++);
    t_float f = atan2f(f1,f2);
    SETFLOAT(m, f);
    m++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}
static void mtx_atan2_float(t_mtx_binmtx *x, t_float f)
{
  t_matrix *m=&x->m, *m2=&x->m2;
  t_atom *ap, *ap2=m2->atombuffer+2;
  int row2, col2, n;

  if (!m2->atombuffer) {
    pd_error(x, "right-hand matrix is missing");
    return;
  }

  row2=atom_getfloat(m2->atombuffer);
  col2=atom_getfloat(m2->atombuffer+1);
  adjustsize(m, row2, col2);
  ap=m->atombuffer+2;

  n=row2*col2;

  while(n--) {
    SETFLOAT(ap, atan2f(f,atom_getfloat(ap2++)));
    ap++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), m->row*m->col+2,
                  m->atombuffer);
}
static void *mtx_atan2_new(t_symbol *s, int argc, t_atom *argv)
{
  if (argc>1) {
    post("mtx_atan2 : extra arguments ignored");
  }
  if (argc) {
    t_mtx_binscalar *x = (t_mtx_binscalar *)pd_new(mtx_atan2scalar_class);
    floatinlet_new(&x->x_obj, &x->f);
    x->f = atom_getfloatarg(0, argc, argv);
    outlet_new(&x->x_obj, 0);
    return(x);
  } else {
    t_mtx_binmtx *x = (t_mtx_binmtx *)pd_new(mtx_atan2_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"), gensym(""));
    outlet_new(&x->x_obj, 0);
    x->m.col = x->m.row =  x->m2.col = x->m2.row = 0;
    x->m.atombuffer = x->m2.atombuffer = 0;
    return(x);
  }
}

void mtx_atan2_setup(void)
{
  mtx_atan2_class = class_new(gensym("mtx_atan2"),
                              (t_newmethod)mtx_atan2_new, (t_method)mtx_binmtx_free,
                              sizeof(t_mtx_binmtx), 0, A_GIMME, 0);
  class_addmethod(mtx_atan2_class, (t_method)mtx_atan2_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addmethod(mtx_atan2_class, (t_method)mtx_bin_matrix2, gensym(""),
                  A_GIMME, 0);
  class_addfloat (mtx_atan2_class, mtx_atan2_float);
  class_addbang  (mtx_atan2_class, mtx_binmtx_bang);

  mtx_atan2scalar_class = class_new(gensym("mtx_atan2"), 0,
                                    (t_method)mtx_binscalar_free,
                                    sizeof(t_mtx_binscalar), 0, 0);
  class_addmethod(mtx_atan2scalar_class, (t_method)mtx_atan2scalar_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_atan2scalar_class, mtx_atan2scalar_list);
  class_addbang  (mtx_atan2scalar_class, mtx_binscalar_bang);



}

void iemtx_atan2_setup(void)
{
  mtx_atan2_setup();
}
