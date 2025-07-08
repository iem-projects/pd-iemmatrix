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

/* mtx_dbtorms: B=dbtorms(A); B[n,m]=dbtorms(A[n,m])  */


void mtx_dbtorms_setup(void)
{
  iemmatrix_unop_setup("mtx_dbtorms", 0, dbtorms, (char*)0);
}

void iemtx_dbtorms_setup(void)
{
  mtx_dbtorms_setup();
}
