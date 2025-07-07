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

/* mtx_resize */

static t_class *mtx_resize_class;
static void mtx_resize_list2(t_matrixobj *x, t_symbol *s, int argc,
                             t_atom *argv)
{
  int r, c;
  (void)s; /* unused */
  if (argc<1) {
    return;
  }
  if (argc>2) {
    pd_error(x,
             "mtx_resize : only rows & cols are needed, skipping the rest");
  }
  if (argc==1) {
    r=c=atom_getfloat(argv++);
  } else {
    r=atom_getfloat(argv++);
    c=atom_getfloat(argv++);
  }

  if (r<0) {
    r=0;
  }
  if (c<0) {
    c=0;
  }

  x->current_row = r;
  x->current_col = c;
}

static void mtx_resize_matrix(t_matrixobj *x, t_symbol *s, int argc,
                              t_atom *argv)
{
  int row=atom_getfloat(argv);
  int col=atom_getfloat(argv+1);
  int r = x->current_row, c = x->current_col;
  int R=0, ROW, COL;
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  if (!r) {
    r=row;
  }
  if (!c) {
    c=col;
  }

  if (r==row && c==col) { /* no need to change */
    outlet_anything(x->x_obj.ob_outlet, gensym("matrix"), argc, argv);
    return;
  }

  x->m.atombuffer=(t_atom *)getbytes((c*r+2)*sizeof(t_atom));
  setdimen(&x->m, r, c);
  matrix_set(&x->m, 0);

  ROW=(r<row)?r:row;
  COL=(c<col)?c:col;
  R=ROW;
  while(R--) {
    memcpy(x->m.atombuffer+2+(ROW-R-1)*c, argv+2+(ROW-R-1)*col,
           COL*sizeof(t_atom));
  }

  matrixobj_bang(x);

  freebytes(x->m.atombuffer, (c*r+2)*sizeof(t_atom));
}

static void *mtx_resize_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_resize_class);
  int c=0, r=0;
  (void)s; /* unused */

  if(argc) {
    if(argc-1) {
      r=atom_getfloat(argv);
      c=atom_getfloat(argv+1);
    } else {
      r=c=atom_getfloat(argv);
    }
    if(c<0) {
      c=0;
    }
    if(r<0) {
      r=0;
    }
  }
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym(""));
  outlet_new(&x->x_obj, 0);
  x->current_row = r;
  x->current_col = c;
  x->m.row = x->m.col= 0;
  x->m.atombuffer  = 0;

  return (x);
}
void mtx_resize_setup(void)
{
  mtx_resize_class = class_new(gensym("mtx_resize"),
                               (t_newmethod)mtx_resize_new,
                               0, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addmethod  (mtx_resize_class, (t_method)mtx_resize_matrix,
                    gensym("matrix"), A_GIMME, 0);
  class_addmethod  (mtx_resize_class, (t_method)mtx_resize_list2,
                    gensym(""), A_GIMME, 0);

}
void iemtx_resize_setup(void)
{
  mtx_resize_setup();
}
