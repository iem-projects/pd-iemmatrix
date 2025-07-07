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

static t_float binop(t_float f1, t_float f2) {
  return f1<=f2;
}


void mtx_le_setup(void)
{
  iemmatrix_binop_setup("mtx_le", 0, binop, "mtx_<=", 0);
}

void iemtx_le_setup(void)
{
  mtx_le_setup();
}
