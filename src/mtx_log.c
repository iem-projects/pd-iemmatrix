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

/* mtx_log: B=log(A); */

static t_float unop(t_float f) {
#if PD_FLOATSIZE == 32
  return logf(f);
#else
  return log(f);
#endif
}


void mtx_log_setup(void)
{
  iemmatrix_unop_setup("mtx_log", unop, 0);
}

void iemtx_log_setup(void)
{
  mtx_log_setup();
}
