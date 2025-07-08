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

/* mtx_not: B=!A; */

#define MTX_ALMOSTZERO 1e-19

static t_float unop(t_float f) {
  return (t_float)(f<MTX_ALMOSTZERO&&f>-MTX_ALMOSTZERO);
}

void mtx_not_setup(void)
{
  iemmatrix_unop_setup("mtx_not", 0, unop, "mtx_!", (char*)0);
}

void iemtx_not_setup(void)
{
  mtx_not_setup();
}
