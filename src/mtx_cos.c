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

/* mtx_cos: B=cos(A); */

static t_float unop(t_float f) {
#if PD_FLOATSIZE == 32
  return cosf(f);
#else
  return cos(f);
#endif
}


void mtx_cos_setup(void)
{
  iemmatrix_unop_setup("mtx_cos", unop, 0);
}

void iemtx_cos_setup(void)
{
  mtx_cos_setup();
}
