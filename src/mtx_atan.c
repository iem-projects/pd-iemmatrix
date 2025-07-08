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

/* mtx_atan: B=atan(A); */

static t_float unop(t_float f) {
#if PD_FLOATSIZE == 32
  return atanf(f);
#else
  return atan(f);
#endif
}


void mtx_atan_setup(void)
{
  iemmatrix_unop_setup("mtx_atan", 0, unop, (char*)0);
}

void iemtx_atan_setup(void)
{
  mtx_atan_setup();
}
