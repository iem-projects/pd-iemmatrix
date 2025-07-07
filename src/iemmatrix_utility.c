/*
 *  iemmatrix_utility
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

/*
  G.Holzmann: this has been in mtx_matrix.c before
              now here should be the shared code !!!
*/
#if defined _WIN32
# include <io.h>
# include <windows.h>
# define snprintf _snprintf
#else
# define _GNU_SOURCE
# include <dlfcn.h>
#endif

#include <unistd.h>

#include "iemmatrix.h"

/* utility functions */

void setdimen(t_matrix *m, int row, int col)
{
  m->col = col;
  m->row = row;
  if(!m->atombuffer) {
    return;
  }
  SETFLOAT(m->atombuffer,   row);
  SETFLOAT(m->atombuffer+1, col);
}

void adjustsize(void*x, t_matrix *m, int desiredRow, int desiredCol)
{
  int col=m->col, row=m->row;

  if (desiredRow<1) {
    pd_error(x, "matrix: cannot make less than 1 rows");
    desiredRow=1;
  }
  if (desiredCol<1) {
    pd_error(x, "matrix: cannot make less than 1 columns");
    desiredCol=1;
  }

  if (col*row!=desiredRow*desiredCol) {
    if(m->atombuffer) {
      freebytes(m->atombuffer, (col*row+2)*sizeof(t_atom));
    }
    m->atombuffer=(t_atom *)getbytes((desiredCol*desiredRow+2)*sizeof(t_atom));
  }

  setdimen(m, desiredRow, desiredCol);
  return;
}

void debugmtx(int argc, t_float *buf, int id)
{
  int i=argc;
  while(i--) {
    int j=argc;
    startpost("debug%d: ", id);
    while(j--) {
      startpost("%f  ", *buf++);
    }
    endpost();
  }
}

t_matrixfloat *matrix2float(t_atom *ap)
{
  int row = atom_getfloat(ap++);
  int col = atom_getfloat(ap++);
  int length = row * col;
  t_matrixfloat *buffer = (t_matrixfloat *)getbytes(sizeof(
                            t_matrixfloat)*length);
  t_matrixfloat *buf = buffer;
  while(length--) {
    *buf++=atom_getfloat(ap++);
  }
  return buffer;
}

void float2matrix(t_atom *ap, t_matrixfloat *buffer)
{
  int row=atom_getfloat(ap++);
  int col=atom_getfloat(ap++);
  int length = row * col;
  t_matrixfloat*buf= buffer;
  while(length--) {
    SETFLOAT(ap, *buf++);
    ap++;
  }
  freebytes(buffer, row*col*sizeof(t_matrixfloat));
}

/* core functions */
void matrix_bang(t_object *x, t_matrix *m)
{
  /* output the matrix */
  if (m->atombuffer) {
    outlet_anything(x->ob_outlet, gensym("matrix"), m->col*m->row+2,
                    m->atombuffer);
  }
}

void matrix_matrix2(void *x, t_matrix *m, int argc, t_atom *argv)
{
  int row, col;
  if(iemmatrix_check(x, 0, argc, argv, 0))return;
  row = atom_getfloat(argv);
  col = atom_getfloat(argv+1);

  /* this is fast and dirty, MAYBE make it slow and clean */
  /* or, to clean matrices, use the mtx_check object */
  if (row*col != m->row*m->col) {
    freebytes(m->atombuffer, m->row*m->col*sizeof(t_atom));
    m->atombuffer = copybytes(argv, (row*col+2)*sizeof(t_atom));
  } else {
    memcpy(m->atombuffer, argv, (row*col+2)*sizeof(t_atom));
  }

  setdimen(m, row, col);
}

void matrixobj_bang(t_matrixobj*x) {
  matrix_bang(&x->x_obj, &x->m);
}
void matrixobj_matrix2(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_matrix2(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_zeros(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_zeros(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_ones(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_ones(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_eye(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_eye(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_egg(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_egg(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_diag(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_diag(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_diegg(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  (void)s;
  matrix_diegg(&x->x_obj, &x->m, argc, argv);
}
void matrixobj_free(t_matrixobj *x)
{
  matrix_free(&x->m);
}
/* set data */

void matrix_set(t_matrix *m, t_float f)
{
  int size = m->col * m->row;
  t_atom *buf=m->atombuffer+2;
  if(m->atombuffer)while(size--) {
      SETFLOAT(&buf[size], f);
    }
}

void matrix_zeros(t_object *x, t_matrix *m, int argc, t_atom *argv)
{
  int col, row;

  switch(argc) {
  case 0: /* zero out the actual matrix */
    matrix_set(m, 0);
    break;
  case 1:
    row=atom_getfloat(argv);
    adjustsize(x, m, row, row);
    matrix_set(m, 0);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
    adjustsize(x, m, row, col);

    matrix_set(m, 0);
  }

  matrix_bang(x, m);
}

void matrix_ones(t_object *x, t_matrix *m, int argc, t_atom *argv)
{
  int col, row;

  switch(argc) {
  case 0: /* zero out the actual matrix */
    matrix_set(m, 1);
    break;
  case 1:
    row=atom_getfloat(argv);
    adjustsize(x, m, row, row);
    matrix_set(m, 1);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
    adjustsize(x, m, row, col);

    matrix_set(m, 1);
  }

  matrix_bang(x, m);
}

void matrix_eye(t_object *x, t_matrix *m, int argc, t_atom *argv)
{
  int col, row;
  int n;

  switch(argc) {
  case 0: /* zero out the actual matrix */
    matrix_set(m, 0);
    break;
  case 1:
    row=atom_getfloat(argv);
    adjustsize(x, m, row, row);
    matrix_set(m, 0);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
    adjustsize(x, m, row, col);
    matrix_set(m, 0);
  }

  col=m->col;
  row=m->row;
  n = (col<row)?col:row;
  while(n--) {
    SETFLOAT(m->atombuffer+2+n*(1+col), 1);
  }

  matrix_bang(x, m);
}

void matrix_egg(t_object *x, t_matrix *m, int argc, t_atom *argv)
{
  int col, row;
  int n;

  switch(argc) {
  case 0: /* zero out the actual matrix */
    matrix_set(m, 0);
    break;
  case 1:
    row=atom_getfloat(argv);
    adjustsize(x, m, row, row);
    matrix_set(m, 0);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
    adjustsize(x, m, row, col);
    matrix_set(m, 0);
  }

  col=m->col;
  row=m->row;
  n = (col<row)?col:row;
  while(n--) {
    SETFLOAT(m->atombuffer+2+(n+1)*(col-1), 1);
  }

  matrix_bang(x, m);
}

void matrix_diag(t_object *x, t_matrix *m, int argc, t_atom *argv)
{
  int col=argc;

  argv+=argc-1;
  if (argc<1) {
    pd_error(x, "matrix: no diagonal present");
    return;
  }
  adjustsize(x, m, argc, argc);
  matrix_set(m, 0);

  while(argc--) {
    SETFLOAT(m->atombuffer+2+argc*(1+col), atom_getfloat(argv--));
  }

  matrix_bang(x, m);
}

void matrix_diegg(t_object *x, t_matrix*m, int argc, t_atom *argv)
{
  int col=argc;

  argv+=argc-1;
  if (argc<1) {
    pd_error(x, "matrix: no dieggonal present");
    return;
  }
  adjustsize(x, m, argc, argc);
  matrix_set(m, 0);

  while(argc--) {
    t_atom *ap=m->atombuffer+2+(argc+1)*(col-1);
    SETFLOAT(ap, atom_getfloat(argv--));
  }

  matrix_bang(x, m);
}


/* destructor */

void matrix_free(t_matrix *m)
{
  freebytes(m->atombuffer, (m->col*m->row+2)*sizeof(t_atom));
  m->atombuffer=0;
  m->col=m->row=0;
}


/* some math */

/*  invert a square matrix (row=col=rowcol) */
/* if "error" is non-NULL, it's content will be set to 0 if the matrix was invertable, else to non-0 */
t_matrixfloat* mtx_doInvert(t_matrixfloat*input, int rowcol, int*err)
{
  /*
   * row==col==rowclo
   * input=t_matrixfloat[row*col]
   * output=t_matrixfloat[row*col]
   */
  int i, k;
  t_matrixfloat *a1, *b1, *a2, *b2;

  int ok=0; /* error counter */

  int col=rowcol, row=rowcol, row2=row*col;
  t_matrixfloat *original=input;
  t_matrixfloat *inverted = 0;

  if(input==0) {
    return 0;
  }

  /* 1a reserve space for the inverted matrix */
  inverted=(t_matrixfloat *)getbytes(sizeof(t_matrixfloat)*row2);
  if(inverted==0) {
    return 0;
  }

  /* 1b make an eye-shaped float-buf for B */
  i=row2;
  b1=inverted;
  while(i--) {
    *b1++=0;
  }
  i=row;
  b1=inverted;
  while(i--) {
    b1[i*(row+1)]=1;
  }

  /* 2. do the Gauss-Jordan */
  for (k=0; k<row; k++) {
    /* adjust current row */
    t_matrixfloat diagel = original[k*(col+1)];
    t_matrixfloat i_diagel = diagel?1./diagel:0;
    if (!diagel) {
      ok++;
    }

    /* normalize current row (set the diagonal-element to 1 */
    a2=original+k*col;
    b2=inverted+k*col;
    i=row;
    while(i--) {
      *a2++*=i_diagel;
      *b2++*=i_diagel;
    }

    /* eliminate the k-th element in each row by adding the weighted normalized row */
    a2=original+k*row;
    b2=inverted+k*row;
    for(i=0; i<row; i++)
      if (i-k) {
        t_matrixfloat f=-*(original+i*row+k);
        int j = row;
        a1=original+i*row;
        b1=inverted+i*row;
        while (j--) {
          *(a1+j)+=f**(a2+j);
          *(b1+j)+=f**(b2+j);
        }
      }
  }
  if(err!=0) {
    *err=ok;
  }

  return inverted;
}

/*  transpose a matrix */
t_matrixfloat*mtx_doTranspose(t_matrixfloat*transposee, int row, int col)
{
  int r,c;
  t_matrixfloat*transposed=0;
  if(!transposee||!row||!col) {
    return 0;
  }
  transposed=(t_matrixfloat*)getbytes(sizeof(t_matrixfloat)*row*col);
  r=row;
  while(r--) {
    c=col;
    while(c--) {
      transposed[c*row+r]=transposee[r*col+c];
    }
  }
  return transposed;
}

/*  multiply matrix A=[rowA*colA] with matrix B=[rowB*colB]; C=A*B; colA=rowB=colArowB */
t_matrixfloat*mtx_doMultiply(int rowA, t_matrixfloat*A, int colArowB,
                             t_matrixfloat*B, int colB)
{
  t_matrixfloat*result=0;
  int r, c, n;

  if(!A || !B || !rowA || !colArowB || !colB) {
    return 0;
  }
  result=(t_matrixfloat*)getbytes(sizeof(t_matrixfloat)*rowA*colB);

  for(r=0; r<rowA; r++) {
    for(c=0; c<colB; c++) {
      t_matrixfloat sum=0.f;
      for(n=0; n<colArowB; n++) {
        sum+=A[colArowB*r+n]*B[colB*n+c];
      }
      result[colB*r+c]=sum;
    }
  }
  return result;
}


void*iemmatrix_getpdfun(const char*name)
{
#ifdef _WIN32
    // get a handle to the module containing the Pd API functions.
    // NB: GetModuleHandle("pd.dll") does not cover all cases.
    HMODULE module;
    if (GetModuleHandleEx(
        GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            (LPCSTR)&pd_typedmess, &module)) {
        return (void *)GetProcAddress(module, name);
    }
#else
    // search recursively, starting from the main program
    return dlsym(dlopen(0, RTLD_NOW), name);
#endif
    return 0;
}

const char*iemmatrix_objname(void*obj) {
  t_object*x = obj;
  t_symbol*s = gensym("");
  if(x && x->te_binbuf) {
    char buf[MAXPDSTRING];
    t_symbol*objsym = atom_getsymbol(binbuf_getvec(x->te_binbuf));
    if(snprintf(buf, MAXPDSTRING, "[%s]: ", objsym->s_name) > 0) {
      buf[MAXPDSTRING-1] = 0;
      s = gensym(buf);
    }
  }
  return s->s_name;
}

const char*iemmatrix_parentabstractionname(const t_glist*ccurrent) {
#if IEMMATRIX_HAVE_G_CANVAS
  t_canvas *canvas = 0;
  t_glist*current = (t_glist*)ccurrent;
  if(0 == current)
    current = canvas_getcurrent();

  canvas = glist_getcanvas(current);
  while(canvas && !canvas_isabstraction(canvas)) {
    canvas=canvas->gl_owner;
  }
  if(canvas && canvas->gl_name) {
    return canvas->gl_name->s_name;
  }
#else
  (void)ccurrent; /* unused */
#endif /* IEMMATRIX_HAVE_G_CANVAS */
  return 0;
}


int iemmatrix_check(void*object, t_symbol*s, int argc, t_atom*argv, unsigned int tests) {
  t_object*x = (t_object*)object;
  const char*objname=iemmatrix_objname(x);

  int row=(argc>1)?atom_getfloat(argv+0):0;
  int col=(argc>1)?atom_getfloat(argv+1):0;
  (void)s; /* unused */

  if (!tests)
    tests =
      IEMMATRIX_CHECK_CRIPPLED
      | IEMMATRIX_CHECK_DIMENSIONS
      | IEMMATRIX_CHECK_SPARSE;

  if ((tests & IEMMATRIX_CHECK_CRIPPLED) && argc<2) {
    pd_error(x, "%scrippled matrix", objname);
    return 1;
  }

  if ((tests & IEMMATRIX_CHECK_DIMENSIONS) && ((col<1)||(row<1))) {
    pd_error(x, "%sinvalid dimensions %dx%d", objname, col, row);
    return 1;
  }

  if ((tests & IEMMATRIX_CHECK_SPARSE)&&(col*row>argc-2)) {
    pd_error(x, "%ssparse matrix not yet supported : use [mtx_check]", objname);
    return 1;
  }
  return 0;
}


void iemmatrix_floats2list(t_atom* dest, const t_float* src, size_t n)
{
  for (; n--; src++, dest++) {
    SETFLOAT (dest, *src);
  }
}
void iemmatrix_floats2list_modulo(t_atom* dest, const t_float* src, size_t n, size_t m)
{
  t_atom *ptr = dest;
  int count1, count2;
  n /= m;
  count1 = m;
  while (count1--)
    for (count2 = n, ptr = dest++; count2--; ptr += m, src++) {
      SETFLOAT(ptr, *src);
    }
}

void iemmatrix_list2floats(t_float* dest, const t_atom* src, size_t n)
{
  while (n--) {
    *dest++ = atom_getfloat (src++);
  }
}

void iemmatrix_list2floats_modulo(t_float* dest, const t_atom* src, size_t n, size_t m)
{
  const t_atom *ptr = src;
  size_t count1, count2;
  n /= m;
  count1 = m;
  while (count1--)
    for (count2 = n, ptr = src++; count2--; ptr += m, dest++) {
      *dest = atom_getfloat (ptr);
    }
}

unsigned int iemmatrix_atom_getuint(const t_atom*a)
{
  t_float f = atom_getfloat(a);
  return (f>0)?(unsigned int)f:0;
}

int iemmatrix_atom_getint(const t_atom*a)
{
  t_float f = atom_getfloat(a);
  return (int)f;
}


struct _iemmatrix_map {
  t_symbol*key;
  void*value;
  struct _iemmatrix_map*next;
};

/* add a new element to the map (pass 'NULL' as <map> to create a new one */
struct _iemmatrix_map*iemmatrix_map_add(struct _iemmatrix_map*map, t_symbol*key, void*value) {
  struct _iemmatrix_map*m = 0;
  for(m=map; m; m=m->next) {
    if(m->key == key) {
      /* m->value = value */
      return map;
    }
  }
  m = getbytes(sizeof(*m));
  m->key = key;
  m->value = value;
  m->next = map;

  return m;
}
void*iemmatrix_map_get(struct _iemmatrix_map*map, t_symbol*key) {
  while(map) {
    if (key == map->key)
      return map->value;
    map = map->next;
  }
  return 0;
}
/* free the map: free_callback() is called for each <value> in the map */
void iemmatrix_map_free(struct _iemmatrix_map*map, void(*free_callback)(void*)) {
  struct _iemmatrix_map*m = 0;
  if(free_callback) {
    for(m = map; m; m=m->next) {
      free_callback(m->value);
    }
  }
  while(map) {
    m = map;
    map = map->next;

    m->key = 0;
    m->value = 0;
    m->next = 0;
    freebytes(m, sizeof(*m));
  }

}
