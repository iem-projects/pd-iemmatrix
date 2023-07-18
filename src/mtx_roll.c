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

/* mtx_roll */
/* roll the rows */
static t_class *mtx_roll_class;

static void mtx_roll_matrix(t_matrix *x, t_symbol *s, int argc,
                            t_atom *argv)
{
  int row, col, colroll;
  t_atom *ap;
  int c;
  if(iemmatrix_check(x, argc, argv, 0))return;
  row=atom_getint(argv++);
  col=atom_getint(argv++);
  colroll = ((int)x->f%col+col)%col;

  adjustsize(x, row, col);
  ap = x->atombuffer+2;

  c=col;
  while(c--) {
    t_atom *in  = argv+col-c-1;
    t_atom *out = ap  +(col-c-1+colroll)%col;
    int r = row;
    while (r--) {
      SETFLOAT(out, atom_getfloat(in));
      out+=col;
      in+=col;
    }

  }

  matrix_bang(x);
}

static void *mtx_roll_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrix *x = (t_matrix *)pd_new(mtx_roll_class);
  floatinlet_new(&x->x_obj, &(x->f));
  outlet_new(&x->x_obj, 0);

  x->f=argc?atom_getfloat(argv):0;
  x->col=x->row=0;
  x->atombuffer=0;
  return (x);
}
void mtx_roll_setup(void)
{
  mtx_roll_class = class_new(gensym("mtx_roll"), (t_newmethod)mtx_roll_new,
                             (t_method)matrix_free, sizeof(t_matrix), 0, A_GIMME, 0);
  class_addbang  (mtx_roll_class, matrix_bang);
  class_addmethod(mtx_roll_class, (t_method)mtx_roll_matrix,
                  gensym("matrix"), A_GIMME, 0);

}
void iemtx_roll_setup(void)
{
  mtx_roll_setup();
}
