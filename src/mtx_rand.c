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

/* mtx_rand */
static t_class *mtx_rand_class;

static void mtx_rand_seed(t_matrixobj *x, t_float f)
{
  x->current_row=f;
}
static int makeseed(void)
{
  static unsigned int random_nextseed = 1489853723;
  random_nextseed = random_nextseed * 435898247 + 938284287;
  return (random_nextseed & 0x7fffffff);
}
static inline t_float getrand(int *val)
{
  t_float f=(((t_float)(((*val=*val*435898247+382842987)&0x7fffffff)
                        -0x40000000))*(t_float)(0.5/0x40000000)+0.5);
  return f;

}
static void mtx_rand_random(t_matrixobj *x)
{
  long size = x->m.row * x->m.col;
  t_atom *ap=x->m.atombuffer+2;
  while(size--) {
    SETFLOAT(ap+size, getrand(&x->current_row));
  }
}

static void mtx_rand_list(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  int row, col;
  (void)s; /* unused */

  switch(argc) {
  case 0:
    return;
  case 1:
    row = col = atom_getfloat(argv);
    break;
  default:
    row = atom_getfloat(argv+0);
    col = atom_getfloat(argv+1);
  }

  adjustsize(x, &x->m, row, col);
  mtx_rand_random(x);
  matrixobj_bang(x);
}
static void mtx_rand_matrix(t_matrixobj *x, t_symbol *s, int argc,
                            t_atom *argv)
{
  (void)s;
  matrix_matrix2(x, &x->m, argc, argv);
  mtx_rand_random(x);
  matrixobj_bang(x);
}
static void mtx_rand_bang(t_matrixobj *x)
{
  if(0==x->m.col || 0==x->m.row) {
    outlet_float(x->x_obj.ob_outlet, getrand(&x->current_row));
  } else {
    mtx_rand_random(x);
    matrixobj_bang(x);
  }
}
static void *mtx_rand_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(mtx_rand_class);
  int row, col;
  (void)s; /* unused */
  outlet_new(&x->x_obj, 0);
  x->m.col=x->m.row=0;
  x->m.atombuffer=0;
  x->current_row=makeseed();

  if (argc) {
    row=atom_getfloat(argv);
    col=(argc>1)?atom_getfloat(argv+1):row;
    adjustsize(x, &x->m, row, col);
    mtx_rand_random(x);
  }
  return (x);
}
void mtx_rand_setup(void)
{
  mtx_rand_class = class_new(gensym("mtx_rand"), (t_newmethod)mtx_rand_new,
                             (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addmethod(mtx_rand_class, (t_method)mtx_rand_matrix,
                  gensym("matrix"), A_GIMME, 0);
  class_addlist  (mtx_rand_class, mtx_rand_list);
  class_addbang  (mtx_rand_class, mtx_rand_bang);

  class_addmethod(mtx_rand_class, (t_method)mtx_rand_seed, gensym("seed"),
                  A_FLOAT, 0);

}

void iemtx_rand_setup(void)
{
  mtx_rand_setup();
}
