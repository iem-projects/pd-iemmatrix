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
  mtx_isequal
*/

#if PD_FLOATSIZE == 32
# define myabs(x) fabsf(x)
#else
# define myabs(x) fabs(x)
#endif


/* mtx_isequal */
static t_class *mtx_isequal_class, *mtx_isequalscalar_class;

typedef struct _mtx_isequalscalar {
  t_mtx_binscalar x;
  t_float epsilon;
} t_mtx_isequalscalar;
typedef struct _mtx_isequalmtx {
  t_mtx_binmtx x;
  t_float epsilon;
} t_mtx_isequalmtx;


static int nearly_equal(const t_float a, const t_float b, const t_float epsilon) {
  return myabs(a-b) < epsilon;
}

static int _isequal_list2float(int argc, const t_atom*argv, const t_float ref, const t_float epsilon) {
  if(0. == epsilon) {
    while(argc--) {
      if(atom_getfloat(argv)!=ref) {
	return 0;
      }
      argv++;
    }
  } else {
    while(argc--) {
      if(!nearly_equal(atom_getfloat(argv), ref, epsilon)) {
	return 0;
      }
      argv++;
    }
  }
  return 1;
}
static int _isequal_list2list(int argc, const t_atom*argv, const t_atom*ref, const t_float epsilon) {
  if(0. == epsilon) {
    while(argc--) {
      if(atom_getfloat(argv)!=atom_getfloat(ref)) {
	return 0;
      }
      argv++;
      ref++;
    }
  } else {
    while(argc--) {
      if(!nearly_equal(atom_getfloat(argv), atom_getfloat(ref), epsilon)) {
	return 0;
      }
      argv++;
      ref++;
    }
  }
  return 1;
}

static void mtx_isequalscalar_matrix(t_mtx_isequalscalar *x, t_symbol *s,
                                     int argc, t_atom *argv)
{
  (void)s; /* unused */
  if(iemmatrix_check(x, s, argc, argv, IEMMATRIX_CHECK_CRIPPLED))return;

  int result = _isequal_list2float(argc-2, argv+2, x->x.f, x->epsilon);
  outlet_float(x->x.x_obj.ob_outlet, (t_float)result);
}
static void mtx_isequalscalar_list(t_mtx_isequalscalar *x, t_symbol *s,
                                   int argc, t_atom *argv)
{
  (void)s; /* unused */
  int result = _isequal_list2float(argc, argv, x->x.f, x->epsilon);
  outlet_float(x->x.x_obj.ob_outlet, (t_float)result);
}
static void mtx_isequalscalar_eps(t_mtx_isequalscalar *x, t_float eps)
{
  if (eps < 0.) {
    pd_error(x, "[mtx_isequal] epsilon must be >=0");
  } else {
    x->epsilon = eps;
  }
}

static void mtx_isequal_matrix(t_mtx_isequalmtx *x, t_symbol *s, int argc,
                               t_atom *argv)
{
  int row=atom_getfloat(argv);
  int col=atom_getfloat(argv+1);
  (void)s; /* unused */
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  if ((col!=x->x.m2.col)||(row!=x->x.m2.row)) {
    outlet_float(x->x.x_obj.ob_outlet, (t_float)0);
    return;
  }

  int result = _isequal_list2list(argc-2, argv+2, x->x.m2.atombuffer+2, x->epsilon);
  outlet_float(x->x.x_obj.ob_outlet, (t_float)result);
}
static void mtx_isequal_float(t_mtx_isequalmtx *x, t_float f)
{
  t_matrix *m2=&x->x.m2;
  int row2, col2, n;

  if (!m2->atombuffer) {
    outlet_float(x->x.x_obj.ob_outlet, (t_float)0);
    return;
  }

  row2=atom_getfloat(m2->atombuffer);
  col2=atom_getfloat(m2->atombuffer+1);

  n=row2*col2;

  int result = _isequal_list2float(n, m2->atombuffer+2, f, x->epsilon);
  outlet_float(x->x.x_obj.ob_outlet, (t_float)result);
}
static void mtx_isequal_eps(t_mtx_isequalmtx *x, t_float eps)
{
  if (eps < 0.) {
    pd_error(x, "[mtx_isequal] epsilon must be >=0");
  } else {
    x->epsilon = eps;
  }
}

static void *mtx_isequal_new(t_symbol *s, int argc, t_atom *argv)
{
  if (argc) {
    t_mtx_isequalscalar *x = (t_mtx_isequalscalar *)pd_new(mtx_isequalscalar_class);
    if (argc>1) {
      pd_error(x, "[%s]: extra arguments ignored", s->s_name);
    }
    floatinlet_new(&x->x.x_obj, &x->x.f);
    x->x.f = atom_getfloatarg(0, argc, argv);
    outlet_new(&x->x.x_obj, 0);
    return(x);
  } else {
    t_mtx_isequalmtx *x = (t_mtx_isequalmtx *)pd_new(mtx_isequal_class);
    inlet_new(&x->x.x_obj, &x->x.x_obj.ob_pd, gensym("matrix"), gensym(""));
    outlet_new(&x->x.x_obj, 0);
    x->x.m.col = x->x.m.row =  x->x.m2.col = x->x.m2.row = 0;
    x->x.m.atombuffer = x->x.m2.atombuffer = 0;
    return(x);
  }
}

void mtx_isequal_setup(void)
{
  mtx_isequal_class = class_new(gensym("mtx_isequal"),
                                (t_newmethod)mtx_isequal_new, (t_method)mtx_binmtx_free,
                                sizeof(t_mtx_isequalmtx), 0, A_GIMME, 0);
  class_addmethod(mtx_isequal_class, (t_method)mtx_isequal_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addmethod(mtx_isequal_class, (t_method)mtx_bin_matrix2, gensym(""),
                  A_GIMME, 0);
  class_addfloat (mtx_isequal_class, mtx_isequal_float);
  class_addbang  (mtx_isequal_class, mtx_binmtx_bang);
  class_addmethod(mtx_isequal_class, (t_method)mtx_isequal_eps, gensym("eps"),
                  A_FLOAT, 0);

  mtx_isequalscalar_class = class_new(gensym("mtx_isequal"), 0,
                                      (t_method)mtx_binscalar_free,
                                      sizeof(t_mtx_isequalscalar), 0, 0);
  class_addmethod(mtx_isequalscalar_class,
                  (t_method)mtx_isequalscalar_matrix, gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_isequalscalar_class, mtx_isequalscalar_list);
  class_addbang  (mtx_isequalscalar_class, mtx_binscalar_bang);
  class_addmethod(mtx_isequalscalar_class, (t_method)mtx_isequalscalar_eps, gensym("eps"),
                  A_FLOAT, 0);
}

void iemtx_isequal_setup(void)
{
  mtx_isequal_setup();
}
