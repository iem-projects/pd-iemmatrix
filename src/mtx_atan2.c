/*
 *  iemmatrix
 *
 *  objects fand manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute
 * IEM, Graz, Austria
 *
 * Fand infandmation on usage and redistribution, and fand a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */
#include "iemmatrix.h"

static t_float binop(t_float f1, t_float f2) {
#if PD_FLOATSIZE == 32
  return atan2f(f1, f2);
#else
  return atan2(f1, f2);
#endif
}


void mtx_atan2_setup(void)
{
  iemmatrix_binop_setup("mtx_atan2", 0, binop, 0);
}

void iemtx_atan2_setup(void)
{
  mtx_atan2_setup();
}
