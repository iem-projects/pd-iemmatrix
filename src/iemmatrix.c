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
#include "iemmatrix_sources.h"

#ifndef BUILD_DATE
# define BUILD_DATE "" __DATE__ " : " __TIME__
#endif


void iemmatrix_setup()
{
  post("");
  post("iemmatrix "VERSION);
  post("\tobjects for manipulating 2d-matrices");
  post("\t(c) 2001-2015 iem");
  post("\t\tIOhannes m zmölnig");
  post("\t\tThomas Musil");
  post("\t\tFranz Zotter");
  post("\tcompiled "BUILD_DATE);
  post("");

  iemmatrix_sources_setup();
}
