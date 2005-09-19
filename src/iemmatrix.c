/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly refering to matlab/octave matrix functions
 *
 * Copyright (c) IOhannes m zm�lnig, forum::f�r::uml�ute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */
#include "iemmatrix.h"

#include "iemmatrix_sources.c"

void iemmatrix_setup(){
  post("");
  post("iemmatrix "VERSION);
  post("\tobjects for manipulating 2d-matrices");
  post("\t(c) IOhannes m zm�lnig, Thomas Musil, Franz Zotter :: iem, 2001-2005");
  post("\tcompiled "__DATE__" : "__TIME__);
  post("");

  iemmatrix_sources_setup();
}
