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

/* mtx_sin: B=sin(A); */

static t_float unop(t_float f) {
#if PD_FLOATSIZE == 32
  return sinf(f);
#else
  return sin(f);
#endif
}


void mtx_sin_setup(void)
{
  iemmatrix_unop_setup("mtx_sin", 0, unop, (char*)0);
}

void iemtx_sin_setup(void)
{
  mtx_sin_setup();
}
