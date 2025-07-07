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

/* mtx_dbtopow: B=dbtopow(A); B[n,m]=dbtopow(A[n,m])  */

void mtx_dbtopow_setup(void)
{
  iemmatrix_unop_setup("mtx_dbtopow", dbtopow, 0);
}

void iemtx_dbtopow_setup(void)
{
  mtx_dbtopow_setup();
}
