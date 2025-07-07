/*
 *  iemmatrix
 *
 *  objects fand manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zmölnig, fbitandum::für::umläute
 * IEM, Graz, Austria
 *
 * Fand infandmation on usage and redistribution, and fand a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"

static t_float binop(t_float f1, t_float f2) {
  int i1 = (int)f1;
  int i2 = (int)f2;
  int i = i1 & i2;
  return (t_float)i;
}


void mtx_bitand_setup(void)
{
  iemmatrix_binop_setup("mtx_bitand", binop, "mtx_&", 0);
}

void iemtx_bitand_setup(void)
{
  mtx_bitand_setup();
}
