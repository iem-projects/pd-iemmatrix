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
#include <stdarg.h>

typedef struct _matrix_binop {
  t_object x_obj;
  t_matrix m;
  t_matrix m1;
  t_matrix m2;
  t_float f2;

  iemmatrix_binopfun_t*fun;

  t_symbol*x_selector;
  int x_argc;
} t_matrix_binop;

typedef struct _binop_ {
  t_class*class;
  t_class*scalarclass;
  iemmatrix_binopfun_t*fun;
} _binop_t;


static struct _iemmatrix_map*s_map = 0;

static void mtx_binop_storematrix(t_matrix_binop*x, t_matrix*m, t_atom*argv) {
  int row = atom_getfloat(argv+0);
  int col = atom_getfloat(argv+1);
  int n = row * col;
  adjustsize(x, m, row, col);
  t_atom *ap = m->atombuffer + 2;
  argv+=2;
  while(n--) {
    t_float f = atom_getfloat(argv++);
    SETFLOAT(ap, f);
    ap++;
  }
}

static void mtx_binop_sendmatrix(t_matrix_binop*x) {
  int argc = x->m.row * x->m.col;
  outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc+2, x->m.atombuffer);
}


/* convert a valid argc/argv matrix to the 'm1' matrix
 * apply a scalar on the matrix
 * if 'reverse' is true, reverse the arguments
 * store the results in x->m
 */
static void mtx_binop_m2f(t_matrix_binop*x, t_matrix*m1, t_float f2, int reverse) {
  iemmatrix_binopfun_t*fun = x->fun;
  t_matrix *m = &x->m; /* output matrix */
  int row = m1->row;
  int col = m1->col;;
  int n = row * col;
  t_atom *ap, *ap1 = m1->atombuffer + 2;
  adjustsize(x, m, row, col);
  ap = m->atombuffer + 2;

  while(n--) {
    t_float f1 = atom_getfloat(ap1++);
    t_float f = reverse?fun(f2, f1):fun(f1, f2);
    SETFLOAT(ap, f);
    ap++;
  }

  mtx_binop_sendmatrix(x);
}
static void mtx_binop_float(t_matrix_binop *x, t_float f1) {
  /* float on left-hand side; apply it to right-hand matrix */
  if(!x->m2.row || !x->m2.col) {
    pd_error(x, "no right-hand matrix to operate on");
    return;
  }
  mtx_binop_m2f(x, &x->m2, f1, 1);
}

static void mtx_binop_bang(t_matrix_binop *x) {
  iemmatrix_binopfun_t*fun = x->fun;
  int row=x->m1.row;
  int col=x->m1.col;
  int n = row * col;

  if (!(x->m2.col && x->m2.row)) {
    mtx_binop_m2f(x, &x->m1, 0., 0);
    return;
  }

  if(x->m2.col==1 && x->m2.row==1) {
    t_float f = atom_getfloat(x->m2.atombuffer+2);
    mtx_binop_m2f(x, &x->m1, f, 0);
    return;
  }

  if(x->m2.row==1 && x->m2.col == x->m1.col) {
    int c, r;
    adjustsize(x, &x->m, row, col);
    t_atom*m = x->m.atombuffer+2;
    t_atom*m1 = x->m1.atombuffer+2;
    for(r=0; r<row; r++) {
      t_atom*m2 = x->m2.atombuffer+2;
      for(c=0; c<col; c++) {
        t_float f = fun(atom_getfloat(m1), atom_getfloat(m2));
        SETFLOAT(m, f);
        m1++;
        m2++;
        m++;
      }
    }
    mtx_binop_sendmatrix(x);
    return;
  }
  if(x->m2.col==1 && x->m2.row == x->m1.row) {
    int c, r;
    adjustsize(x, &x->m, row, col);
    t_atom*m = x->m.atombuffer+2;
    t_atom*m1 = x->m1.atombuffer+2;
    t_atom*m2 = x->m2.atombuffer+2;
    for(r=0; r<row; r++) {
      t_float f2 = atom_getfloat(m2);
      for(c=0; c<col; c++) {
        t_float f = fun(atom_getfloat(m1), f2);
        SETFLOAT(m, f);
        m++;
        m1++;
      }
      m2++;
    }

    mtx_binop_sendmatrix(x);
    return;
  }

  if((x->m1.col != x->m2.col) || (x->m1.row != x->m2.row)) {
    pd_error(x, "matrix dimensions do not match");
    return;
  }
  adjustsize(x, &x->m, row, col);
  t_atom *m  = x->m .atombuffer+2;
  t_atom *m1 = x->m1.atombuffer+2;
  t_atom *m2 = x->m2.atombuffer+2;

  while(n--) {
    t_float f = fun(atom_getfloat(m1++), atom_getfloat(m2++));
    SETFLOAT(m, f);
    m++;
  }

  mtx_binop_sendmatrix(x);

}

static void mtx_binop_matrix(t_matrix_binop *x,
                             t_symbol *s, int argc, t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  x->x_selector = s;
  x->x_argc = argc;
  mtx_binop_storematrix(x, &x->m1, argv);
  mtx_binop_bang(x);
}

static void mtx_binop_matrix2(t_matrix_binop *x,
                              t_symbol *s, int argc, t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  mtx_binop_storematrix(x, &x->m2, argv);
}

static void mtx_binopscalar_bang(t_matrix_binop *x) {
  iemmatrix_binopfun_t*fun = x->fun;
  int row=x->m1.row;
  int col=x->m1.col;
  int n = row * col;

  if(!x->x_argc) {
    outlet_bang(x->x_obj.ob_outlet);
    return;
  }

  adjustsize(x, &x->m, row, col);
  t_atom *m  = x->m .atombuffer+2;
  t_atom *m1 = x->m1.atombuffer+2;
  t_float f2 = x->f2;

  while(n--) {
    t_float f = fun(atom_getfloat(m1++), f2);
    SETFLOAT(m, f);
    m++;
  }

  if(x->x_selector == gensym("matrix"))
    mtx_binop_sendmatrix(x);
  else {
    outlet_anything(x->x_obj.ob_outlet, x->x_selector, x->x_argc, x->m.atombuffer+2);
  }
}

static void mtx_binopscalar_matrix(t_matrix_binop *x,
                                   t_symbol *s, int argc, t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  x->x_selector = s;
  x->x_argc = argc;
  mtx_binop_storematrix(x, &x->m1, argv);
  mtx_binopscalar_bang(x);
}

static void mtx_binopscalar_list(t_matrix_binop *x, t_symbol *s, int argc,
                                 t_atom *argv)
{
  if(s)
    x->x_selector = s;
  else if (1 == argc)
    x->x_selector = gensym("float");
  else
    x->x_selector = gensym("list");
  x->x_argc = argc;

  adjustsize(x, &x->m1, 1, argc);
  t_atom *m1 = x->m1.atombuffer + 2;
  while(argc--) {
    t_float f = atom_getfloat(argv++);
    SETFLOAT(m1, f);
    m1++;
  }
  mtx_binopscalar_bang(x);
}

static void *mtx_binop_new(t_symbol*s, int argc, t_atom*argv)
{
  t_matrix_binop *x = 0;
  /* element cos */
  _binop_t *binop=iemmatrix_map_get(s_map, s);
  if(!binop)
    return x;
  t_class*cls = binop->class;
  t_class*scalarcls = binop->scalarclass;
  iemmatrix_binopfun_t *fun = binop->fun;
  if(!cls && !scalarcls) {
    pd_error(x, "no class for '%s'", s->s_name);
    return x;
  }
  if(!fun) {
    pd_error(x, "no function for '%s'", s->s_name);
    return x;
  }

  if(argc) {
    x = (t_matrix_binop *)pd_new(scalarcls);
    floatinlet_new(&x->x_obj, &x->f2);
    x->f2 = atom_getfloatarg(0, argc, argv);
  } else {
    x = (t_matrix_binop *)pd_new(cls);
    inlet_new(&x->x_obj, &x->x_obj.te_g.g_pd, gensym("matrix"), gensym(""));
    x->m1.col = x->m1.row = x->m2.col = x->m2.row = 0;
    x->m1.atombuffer = x->m2.atombuffer = 0;
  }

  x->fun = fun;

  outlet_new(&x->x_obj, 0);
  if (argc>1) {
    pd_error(x, "%s : extra arguments ignored", (s && s->s_name) ? s->s_name : "matrix object");
  }
  return(x);
}
static void mtx_binop_free(t_matrix_binop*x) {
  matrix_free(&x->m);
  matrix_free(&x->m1);
  matrix_free(&x->m2);
}


void iemmatrix_binop_setup(const char*classname, const char*helpname, iemmatrix_binopfun_t*fun, ...) {
  if(!classname || !fun)
    return;
  t_class*cls, *scalarcls;
  t_symbol*s = gensym(classname);
  t_symbol*help = helpname?gensym(helpname):0;
  va_list ap;
  if(s == help)
    help = 0;

  cls = class_new(s,
                  (t_newmethod)mtx_binop_new,
                  (t_method)mtx_binop_free,
                  sizeof(t_matrix_binop),
                  0,
                  A_GIMME, A_NULL);
  class_addmethod(cls, (t_method)mtx_binop_matrix, gensym("matrix"), A_GIMME, A_NULL);
  class_addmethod(cls, (t_method)mtx_binop_matrix2, gensym(""), A_GIMME, A_NULL);
  class_addfloat(cls, (t_method)mtx_binop_float);
  class_addbang(cls, (t_method)mtx_binop_bang);

  if(help)
    class_sethelpsymbol(cls, help);

  //class_addlist(cls, (t_method)mtx_binop_list);

  scalarcls = class_new(s,
                        0, /* no constructor */
                        (t_method)mtx_binop_free,
                        sizeof(t_matrix_binop),
                        0, A_NULL);
  class_addmethod(scalarcls, (t_method)mtx_binopscalar_matrix, gensym("matrix"), A_GIMME, A_NULL);
  class_addlist(scalarcls, (t_method)mtx_binopscalar_list);
  class_addbang(scalarcls, (t_method)mtx_binopscalar_bang);
  if(help)
    class_sethelpsymbol(scalarcls, help);

  _binop_t*binop = (_binop_t*)getbytes(sizeof(*binop));
  binop->class = cls;
  binop->scalarclass = scalarcls;
  binop->fun = fun;

  s_map = iemmatrix_map_add(s_map, s, binop);

  va_start(ap, fun);
  while(1) {
    const char*alias = va_arg(ap, char*);
    if(!alias)
      break;
    s = gensym(alias);
    class_addcreator((t_newmethod)mtx_binop_new, s, A_GIMME, A_NULL);
    s_map = iemmatrix_map_add(s_map, s, binop);
  }
  va_end(ap);
}
