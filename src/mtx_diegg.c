/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly refering to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zm�lnig, forum::f�r::uml�ute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */
#include "iemmatrix.h"
static t_class *mtx_diegg_class;
static void mtx_diegg_matrix(t_matrix *x, t_symbol *s, int argc, t_atom *argv)
{
  int row=atom_getfloat(argv++);
  int col=atom_getfloat(argv++);
  int length=(col<row)?col:row, n=length;
  t_atom *ap = (t_atom *)getbytes(length * sizeof(t_atom)), *dummy=ap;
  if(row*col>argc-2)post("mtx_diegg: sparse matrices not yet supported : use \"mtx_check\"");
  else {
    for(n=0;n<length;n++, dummy++)SETFLOAT(dummy, atom_getfloat(argv+(n-1)*(col-1)));
    outlet_list(x->x_obj.ob_outlet, gensym("diegg"), length, ap);
  }
  freebytes(ap, (length * sizeof(t_atom)));
}
static void *mtx_diegg_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrix *x = (t_matrix *)pd_new(mtx_diegg_class);
  outlet_new(&x->x_obj, 0);
  x->row = x->col = 0;
  x->atombuffer   = 0;

  if(!argc)return(x);
  x->atombuffer = (t_atom *)getbytes((argc*argc+2)*sizeof(t_atom));
  setdimen(x, argc, argc);
  matrix_set(x, 0);
  argv+=argc-1;
  while(argc--)SETFLOAT(x->atombuffer+2+argc*(1+x->col), atom_getfloat(argv--));

  return (x);
}
void mtx_diegg_setup(void)
{
  mtx_diegg_class = class_new(gensym("mtx_diegg"), (t_newmethod)mtx_diegg_new, 
			     (t_method)matrix_free, sizeof(t_matrix), 0, A_GIMME, 0);
  class_addlist  (mtx_diegg_class, matrix_diegg);
  class_addbang  (mtx_diegg_class, matrix_bang);
  class_addmethod(mtx_diegg_class, (t_method)mtx_diegg_matrix, gensym("matrix"), A_GIMME, 0);

}
void iemtx_diegg_setup(void){
  mtx_diegg_setup();
}
