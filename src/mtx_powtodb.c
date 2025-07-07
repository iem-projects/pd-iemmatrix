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
/* mtx_powtodb: B=powtodb(A); B[n,m]=powtodb(A[n,m])  */


void mtx_powtodb_setup(void)
{
  iemmatrix_unop_setup("mtx_powtodb", 0, powtodb, 0);
}

void iemtx_powtodb_setup(void)
{
  mtx_powtodb_setup();
}
