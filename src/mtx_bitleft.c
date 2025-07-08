/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zmölnig, fbitleftum::für::umläute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"

static t_float binop(t_float f1, t_float f2) {
  int i1 = (int)f1;
  int i2 = (int)f2;
  int i = i1 << i2;
  return (t_float)i;
}


void mtx_bitleft_setup(void)
{
  iemmatrix_binop_setup("mtx_bitleft", 0, binop, "mtx_<<", (char*)0);
}

void iemtx_bitleft_setup(void)
{
  mtx_bitleft_setup();
}
