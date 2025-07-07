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

typedef struct _matrix_unop {
  t_object x_obj;
  t_matrix m;
  iemmatrix_unopfun_t*fun;
  t_symbol*x_selector;
  int x_argc;
} t_matrix_unop;

typedef struct _unop_ {
  t_class*class;
  iemmatrix_unopfun_t*fun;
} _unop_t;


static struct _iemmatrix_map*s_map;

static void mtx_unop_bang(t_matrix_unop *x) {
  if(!x->x_argc)
    outlet_bang(x->x_obj.ob_outlet);
  else
    outlet_anything(x->x_obj.ob_outlet, x->x_selector, x->x_argc, x->m.atombuffer);
}


static void mtx_unop_matrix(t_matrix_unop *x,
                            t_symbol *s, int argc, t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  int row=atom_getint(argv++);
  int col=atom_getint(argv++);
  int n = argc-2;


  adjustsize(x, &x->m, row, col);
  t_atom *m =  x->m.atombuffer+2;
  iemmatrix_unopfun_t*fun = x->fun;

  while(n--) {
    t_float f = fun(atom_getfloat(argv++));
    SETFLOAT(m, f);
    m++;
  }

  x->x_selector = s;
  x->x_argc = argc;
  mtx_unop_bang(x);
}

static void mtx_unop_list(t_matrix_unop *x, t_symbol *s, int argc,
                          t_atom *argv)
{
  int n=argc;
  t_atom *m;
  (void)s; /* unused */
  iemmatrix_unopfun_t*fun = x->fun;

  adjustsize(x, &x->m, 1, argc);
  m = x->m.atombuffer;

  while(n--) {
    t_float f = fun(atom_getfloat(argv++));
    SETFLOAT(m, f);
    m++;
  }

  if(s)
    x->x_selector = s;
  else if (1 == argc)
    x->x_selector = gensym("float");
  else
    x->x_selector = gensym("list");
  x->x_argc = argc;
  mtx_unop_bang(x);
}

static void *mtx_unop_new(t_symbol*s, int argc, t_atom*argv)
{
  /* element cos */
  _unop_t *unop=iemmatrix_map_get(s_map, s);
  if(!unop)
    return 0;
  (void)argc;(void)argv; /* unused */
  t_class*cls = unop->class;
  iemmatrix_unopfun_t *fun = unop->fun;
  if(!cls) {
    pd_error(0, "no class for '%s'", s->s_name);
    return 0;
  }
  if(!fun) {
    pd_error(0, "no function for '%s'", s->s_name);
    return 0;
  }

  t_matrix_unop *x = (t_matrix_unop *)pd_new(cls);
  x->fun = fun;
  outlet_new(&x->x_obj, 0);
  return(x);
}


void iemmatrix_unop_setup(const char*classname, const char*helpname, iemmatrix_unopfun_t*fun, ...) {
  if(!classname || !fun)
    return;
  t_class*cls;
  t_symbol*s = gensym(classname);
  va_list ap;

  _unop_t*unop = (_unop_t*)getbytes(sizeof(*unop));

  cls = class_new(s,
                  (t_newmethod)mtx_unop_new,
                  0,
                  sizeof(t_matrix_unop),
                  0,
                  A_GIMME, A_NULL);

  class_addbang(cls, (t_method)mtx_unop_bang);
  class_addlist(cls, (t_method)mtx_unop_list);
  class_addmethod(cls, (t_method)mtx_unop_matrix, gensym("matrix"), A_GIMME, A_NULL);

  if(helpname) {
    t_symbol*help = gensym(helpname);
    if (help != s)
      class_sethelpsymbol(cls, help);
  }

  unop->class = cls;
  unop->fun = fun;

  s_map = iemmatrix_map_add(s_map, s, unop);

  va_start(ap, fun);
  while(1) {
    const char*alias = va_arg(ap, char*);
    if(!alias)
      break;
    s = gensym(alias);
    class_addcreator((t_newmethod)mtx_unop_new, s, A_GIMME, A_NULL);
    s_map = iemmatrix_map_add(s_map, s, unop);
  }
  va_end(ap);
}
