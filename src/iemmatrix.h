/* ************************************* */
/* iemmatrix                             */
/* ************************************* */
/*  objects for simple matrix operations */
/* ************************************* */

/*
 * IEMMATRIX is a runtime-library
 * for miller s. puckette's realtime-computermusic-software "pure data"
 * therefore you NEED "pure data" to make any use of the IEMMATRIX external
 * (except if you want to use the code for other things)
 *
 * you can get "pure data" at
 *   http://pd.iem.at
 *   ftp://iem.at/pd
 */

/*
 * Copyright (c) Thomas Musil; IEM KUG, Graz, Austria; 2001-2005
 * Copyright (c) IOhannes m zmölnig (forum::für::umläute), IEM KUG, Graz, Austria; 2001-2005
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 *
 * "pure data" has it's own license, that comes shipped with "pure data".
 *
 * there are ABSOLUTELY NO WARRANTIES for anything
 */

#ifndef INCLUDE_IEMMATRIX_H__
#define INCLUDE_IEMMATRIX_H__


#ifdef __WIN32__
/* MinGW automatically defines __WIN32__
 * other windows compilers might have to define it by hand
 */

# pragma warning( disable : 4244 )
# pragma warning( disable : 4305 )

#endif /* __WIN32__ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif /* HAVE_CONFIG_H */


#ifdef HAVE_M_PD_H
# include <m_pd.h>
#elif defined HAVE_PD_M_PD_H
# include <pd/m_pd.h>
#else
# include <m_pd.h>
#endif


#ifdef HAVE_G_CANVAS_H
# include <g_canvas.h>
# define IEMMATRIX_HAVE_G_CANVAS 1
#elif defined HAVE_PD_G_CANVAS_H
# include <pd/g_canvas.h>
# define IEMMATRIX_HAVE_G_CANVAS 1
#else
# define IEMMATRIX_HAVE_G_CANVAS 0
#endif

#ifndef VERSION
#ifdef PACKAGE_VERSION
# define VERSION PACKAGE_VERSION
#else
# define VERSION "(unknown)"
#endif
#endif /* VERSION */

#include <math.h>
#include <string.h>

#ifdef __WIN32__
#ifndef fabsf
# define fabsf fabs
#endif
#ifndef sqrtf
# define sqrtf sqrt
#endif
#ifndef powf
# define powf pow
#endif
#ifndef atanf
# define atanf atan
#endif
#endif

#ifdef __APPLE__
# include <AvailabilityMacros.h>
# if defined (MAC_OS_X_VERSION_10_3) && MAC_OS_X_VERSION_MAX_ALLOWED  >= MAC_OS_X_VERSION_10_3
# else
/* float intrinsics not in math.h, so we define them here */
#  define sqrtf(v)    (float)sqrt((double)(v))
#  define cosf(v)     (float)cos((double)(v))
#  define sinf(v)     (float)sin((double)(v))
#  define tanf(v)     (float)tan((double)(v))
#  define logf(v)     (float)log((double)(v))
#  define expf(v)     (float)exp((double)(v))
#  define atan2f(v,p) (float)atan2((double)(v), (double)(p))
#  define powf(v,p)   (float)pow((double)(v), (double)(p))
# endif
#endif

typedef double t_matrixfloat;

/* the main class...*/
typedef struct _matrix {
  t_object x_obj;

  int      row;
  int      col;

  t_atom *atombuffer;

  int     current_row,
          current_col;  /* this makes things easy for the mtx_row & mtx_col...*/
  t_float f;

  t_canvas *x_canvas; /* needed for file-reading */
  t_outlet *x_outlet; /* just in case somebody wants an outlet */
} t_matrix;

typedef struct _mtx_binscalar {
  t_object x_obj;

  t_matrix m; /* the output matrix */
  t_float f;  /* the second input */
} t_mtx_binscalar;

typedef struct _mtx_binmtx {
  t_object x_obj;

  t_matrix m;  /* the output matrix */
  t_matrix m2; /* the second input */
} t_mtx_binmtx;


/*
  G.Holzmann: the following is now implemented
              in iemmatrix_utility.c
*/

void matrix_free(t_matrix*x);

/* utility function */
void setdimen(t_matrix *x, int row, int col);
void adjustsize(t_matrix *x, int desiredRow, int desiredCol);
void debugmtx(int argc, t_float *buf, int id);
t_matrixfloat *matrix2float(t_atom *ap);
void float2matrix(t_atom *ap, t_matrixfloat *buffer);
/* basic matrix message checks: returns 1 if matrix message is invalid, 0 otherwise ; prints error messages
 * tests specifies the tests to run (0==all tests)
 */
enum {
  IEMMATRIX_CHECK_CRIPPLED   = 1<<0,
  IEMMATRIX_CHECK_DIMENSIONS = 1<<1,
  IEMMATRIX_CHECK_SPARSE     = 1<<2,
  IEMMATRIX_CHECK_ALL        = 0
};
int iemmatrix_check(void*object, int argc, t_atom*argv, unsigned int tests);

/* get a (decorated) objname: returns either "[<foo>]: " or ""
 * the returned string MUST NOT be freed by the caller
 * 'object' must be of type 't_object*' (or derived)
 */
const char*iemmatrix_objname(void*object);

/* get the name of the abstraction containing the object
 * you can use current=0 when calling this in the constructor
 */
const char*iemmatrix_parentabstractionname(const t_glist*current);

/* get a Pd function by name */
void*iemmatrix_getpdfun(const char*name);


/* basic I/O functions */
void matrix_bang(t_matrix *x); /* output the matrix stored in atombuffer */
void matrix_matrix2(t_matrix *x, t_symbol *s, int argc,
                    t_atom *argv); /* store the matrix in atombuffer */

/* set data */
void matrix_set(t_matrix *x, t_float f); /* set the entire matrix to "f" */
void matrix_zeros(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_ones(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_eye(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_egg(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_diag(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_diegg(t_matrix *x, t_symbol *s, int argc, t_atom *argv);

/* get/set data */
void matrix_row(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_col(t_matrix *x, t_symbol *s, int argc, t_atom *argv);
void matrix_element(t_matrix *x, t_symbol *s, int argc, t_atom *argv);


/*
  in iemmatrix_binops.c
*/

void mtx_bin_matrix2(t_mtx_binmtx *x, t_symbol *s, int argc, t_atom *argv);
void mtx_binmtx_bang(t_mtx_binmtx *x);
void mtx_binmtx_free(t_mtx_binmtx *x);
void mtx_binscalar_bang(t_mtx_binscalar *x);
void mtx_binscalar_free(t_mtx_binscalar *x);


/* some math */

/*  invert a square matrix (row=col=rowcol) */
/* if "error" is non-NULL, it's content will be set to 0 if the matrix was invertable, else to non-0 */
t_matrixfloat*mtx_doInvert(t_matrixfloat*input, int rowcol, int*error);
/*  transpose a matrix */
t_matrixfloat*mtx_doTranspose(t_matrixfloat*output, int row, int col);
/*  multiply matrix A=[rowA*colA] with matrix B=[rowB*colB]; C=A*B; colA=rowB=colArowB */
t_matrixfloat*mtx_doMultiply(int rowA, t_matrixfloat*A, int colArowB,
                             t_matrixfloat*B, int colB);
/* version comparison */
#define KERNEL_VERSION(a, b, c)  ((a * 1000) + b) * 1000 + c

/* for debugging purposes */
#define MARK	startpost("MARK[%s:%d@%s]", __FILE__, __LINE__, __FUNCTION__), post


#endif
