/* ************************************* */
/* iemmatrix                             */
/* ************************************* */
/*  objects for simple matrix operations */
/* ************************************* */

/*
 * Copyright (c) IOhannes m zmölnig (forum::für::umläute), IEM KUG, Graz, Austria; 2025
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 * there are ABSOLUTELY NO WARRANTIES for anything
 */
#ifndef _iemmatrix_stub_h_
#define _iemmatrix_stub_h_
#include <stddef.h>

typedef double (*t_besselfn)(int n, double x);

typedef void *(*t_stub_alloc)(const size_t);
typedef void (*t_stub_free)(void *);
typedef void *(*t_stub_alloc2)(const size_t, const size_t);
typedef void *(*t_stub_calloc)(const size_t, const size_t);
typedef double (*t_stub_get)(void *, const size_t, const size_t);
typedef void *(*t_stub_set)(void *, const size_t, const size_t, const double);
#define IEMMATRIX_DECLARE_ALLOCFREE_STUB(var) \
  static t_stub_alloc var##_alloc = 0;        \
  static t_stub_free var##_free = 0;
#define IEMMATRIX_DECLARE_ALLOCFREE2_STUB(var) \
  static t_stub_alloc2 var##_alloc = 0;        \
  static t_stub_calloc var##_calloc = 0;       \
  static t_stub_get var##_get = 0;             \
  static t_stub_set var##_set = 0;             \
  static t_stub_free var##_free = 0;

struct _class;
void *iemmatrix_get_stub(const char *name, struct _class *);

#endif /*  _iemmatrix_stub_h_ */
