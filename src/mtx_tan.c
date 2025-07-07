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

/* mtx_tan: B=tan(A); */

static t_float unop(t_float f) {
#if PD_FLOATSIZE == 32
  return tanf(f);
#else
  return tan(f);
#endif
}


void mtx_tan_setup(void)
{
  iemmatrix_unop_setup("mtx_tan", unop, 0);
}

void iemtx_tan_setup(void)
{
  mtx_tan_setup();
}
