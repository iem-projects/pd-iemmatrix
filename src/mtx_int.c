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

/* mtx_int: B=int(A); */


static t_float unop(t_float f) {
  int i = (int)f;
  return (t_float)i;
}


void mtx_int_setup(void)
{
  iemmatrix_unop_setup("mtx_int", 0, unop, 0);
}

void iemtx_int_setup(void)
{
  mtx_int_setup();
}
