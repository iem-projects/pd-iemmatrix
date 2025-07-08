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

/*
  mtx_mul
  mtx_*
  mtx_.*
  mtx_./

  matrix multiplication

*/


/* mtx_mul */
static t_class *mtx_mul_class, *mtx_mulscalar_class;

static void mtx_mul_matrix(t_mtx_binmtx *x, t_symbol *s, int argc,
                           t_atom *argv)
{
  t_matrix *m=&x->m, *m2=&x->m2;
  t_atom *ap, *ap1=argv+2, *ap2=m2->atombuffer+2;
  int row=atom_getfloat(argv), col=atom_getfloat(argv+1);
  int row2, col2, n, r, c;

  if (!m2->atombuffer) {
    pd_error(x, "[mtx_*]: right-hand matrix is missing");
    return;
  }
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  row2=atom_getfloat(m2->atombuffer);
  col2=atom_getfloat(m2->atombuffer+1);

  if (col!=row2) {
    pd_error(x, "[mtx_*]: matrix dimensions do not match !");
    return;
  }

  adjustsize(x, m, row, col2);
  ap=m->atombuffer+2;

  for(r=0; r<row; r++)
    for(c=0; c<col2; c++) {
      t_matrixfloat sum = 0.f;
      for(n=0; n<col; n++) {
        sum+=(t_matrixfloat)atom_getfloat(ap1+col*r+n)*atom_getfloat(
               ap2+col2*n+c);
      }
      SETFLOAT(ap+col2*r+c,sum);
    }
  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), m->row*m->col+2,
                  m->atombuffer);
}

static void mtx_mul_float(t_mtx_binmtx *x, t_float f)
{
  t_matrix *m=&x->m, *m2=&x->m2;
  t_atom *ap, *ap2=m2->atombuffer+2;
  int row2, col2, n;

  if (!m2->atombuffer) {
    pd_error(x, "[mtx_*]: right-hand matrix is missing");
    return;
  }

  row2=atom_getfloat(m2->atombuffer);
  col2=atom_getfloat(m2->atombuffer+1);
  adjustsize(x, m, row2, col2);
  ap=m->atombuffer+2;

  n=row2*col2;

  while(n--) {
    SETFLOAT(ap, f*atom_getfloat(ap2++));
    ap++;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), m->row*m->col+2,
                  m->atombuffer);
}

static void mtx_mulscalar_matrix(t_mtx_binscalar *x, t_symbol *s, int argc,
                                 t_atom *argv)
{
  int row, col;
  int n=argc-2;
  t_atom *m;
  t_float factor = x->f;

  if(iemmatrix_check(x, s, argc, argv, IEMMATRIX_CHECK_CRIPPLED))return;
  row=atom_getfloat(argv++);
  col=atom_getfloat(argv++);

  adjustsize(x, &x->m, row, col);
  m = x->m.atombuffer+2;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = atom_getfloat(argv++)*factor;
  }

  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc,
                  x->m.atombuffer);
}
static void mtx_mulscalar_list(t_mtx_binscalar *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  int n=argc;
  t_atom *m;
  t_float factor = x->f;
  (void)s; /* unused */
  adjustsize(x, &x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    m->a_type = A_FLOAT;
    (m++)->a_w.w_float = atom_getfloat(argv++)*factor;
  }
  outlet_list(x->x_obj.ob_outlet, gensym("list"), argc, x->m.atombuffer);
}

static void *mtx_mul_new(t_symbol *s, int argc, t_atom *argv)
{
  if (argc>1) {
    pd_error(0, "[%s]: extra arguments ignored", s->s_name);
  }
  if (argc) {
    t_mtx_binscalar *x = (t_mtx_binscalar *)pd_new(mtx_mulscalar_class);
    floatinlet_new(&x->x_obj, &x->f);
    x->f = atom_getfloatarg(0, argc, argv);
    outlet_new(&x->x_obj, 0);
    return(x);
  } else {
    t_mtx_binmtx *x = (t_mtx_binmtx *)pd_new(mtx_mul_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"), gensym(""));
    outlet_new(&x->x_obj, 0);
    x->m.col = x->m.row = x->m2.col = x->m2.row = 0;
    x->m.atombuffer = x->m2.atombuffer = 0;
    return (x);
  }
}

static t_float binop(t_float f1, t_float f2) {
  return f1 * f2;
}

void mtx_mul_setup(void)
{
  mtx_mul_class = class_new(gensym("mtx_mul"), (t_newmethod)mtx_mul_new,
                            (t_method)mtx_binmtx_free,
                            sizeof(t_mtx_binmtx), 0, A_GIMME, 0);
  class_addcreator((t_newmethod)mtx_mul_new, gensym("mtx_*"), A_GIMME,0);
  class_addmethod(mtx_mul_class, (t_method)mtx_mul_matrix, gensym("matrix"),
                  A_GIMME, 0);
  class_addmethod(mtx_mul_class, (t_method)mtx_bin_matrix2, gensym(""),
                  A_GIMME, 0);
  class_addfloat (mtx_mul_class, mtx_mul_float);
  class_addbang  (mtx_mul_class, mtx_binmtx_bang);

  mtx_mulscalar_class = class_new(gensym("mtx_mul"), 0,
                                  (t_method)mtx_binscalar_free,
                                  sizeof(t_mtx_binscalar), 0, 0);

  class_addmethod(mtx_mulscalar_class, (t_method)mtx_mulscalar_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_mulscalar_class, mtx_mulscalar_list);
  class_addbang  (mtx_mulscalar_class, mtx_binscalar_bang);

  iemmatrix_binop_setup("mtx_.*", "mtx_mul", binop, (char*)0);
}

void iemtx_mul_setup(void)
{
  mtx_mul_setup();
}
