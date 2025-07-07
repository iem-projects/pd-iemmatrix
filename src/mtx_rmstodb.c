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
/* mtx_rmstodb: B=rmstodb(A); B[n,m]=rmstodb(A[n,m])  */


void mtx_rmstodb_setup(void)
{
  iemmatrix_unop_setup("mtx_rmstodb", 0, rmstodb, 0);
}

void iemtx_rmstodb_setup(void)
{
  mtx_rmstodb_setup();
}
