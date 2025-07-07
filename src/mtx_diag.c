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

/* mtx_diag */
static t_class *mtx_diag_class;
static void mtx_diag_matrix(t_matrixobj *x, t_symbol *s, int argc,
                            t_atom *argv)
{
  int row, col, length, n;
  t_atom *ap = 0, *dummy=0;
  (void)s; /* unused */
  if(iemmatrix_check(x, s, argc, argv, 0))return;

  row=atom_getfloat(argv++);
  col=atom_getfloat(argv++);
  length=(col<row)?col:row;
  n=length;

  ap=(t_atom *)getbytes(length * sizeof(t_atom));
  dummy=ap;
  for(n=0; n<length; n++, dummy++) {
    SETFLOAT(dummy, atom_getfloat(argv+n*(col+1)));
  }
  outlet_list(x->x_obj.ob_outlet, gensym("diag"), length, ap);
  freebytes(ap, (length * sizeof(t_atom)));
}
static void *mtx_diag_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_diag_class);
  (void)s; /* unused */
  outlet_new(&x->x_obj, 0);
  x->m.row = x->m.col = 0;
  x->m.atombuffer   = 0;

  if(!argc) {
    return(x);
  }
  x->m.atombuffer = (t_atom *)getbytes((argc*argc+2)*sizeof(t_atom));
  setdimen(&x->m, argc, argc);
  matrix_set(&x->m, 0);
  argv+=argc-1;
  while(argc--) {
    SETFLOAT(x->m.atombuffer+2+argc*(1+x->m.col), atom_getfloat(argv--));
  }

  return (x);
}
void mtx_diag_setup(void)
{
  mtx_diag_class = class_new(gensym("mtx_diag"), (t_newmethod)mtx_diag_new,
                             (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addlist  (mtx_diag_class, matrixobj_diag);
  class_addbang  (mtx_diag_class, matrixobj_bang);
  class_addmethod(mtx_diag_class, (t_method)mtx_diag_matrix,
                  gensym("matrix"), A_GIMME, 0);

}
void iemtx_diag_setup(void)
{
  mtx_diag_setup();
}
