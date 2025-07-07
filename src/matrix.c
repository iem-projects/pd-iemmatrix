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

/*
  matrix : basic object : create and store matrices
  mtx    : alias for matrix
*/

#include "iemmatrix.h"
#include <stdio.h>
#ifdef _WIN32
/* or should we use the security enhanced _snprintf_s() ?? */
# define snprintf _snprintf
#endif


/* -------------------- matrix ------------------------------ */
/*
  G.Holzmann: see iemmatrix_utility.c for the functions
              you don't find here ... ;)
*/

static t_class *matrix_class;

static void matrix_matrix(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  if(iemmatrix_check(x, s, argc, argv, 0))return;
  matrix_matrix2(x, &x->m, argc, argv);
  matrixobj_bang(x);
}

static void matrix_float(t_matrixobj *x, t_float f)
{
  matrix_set(&x->m, f);
  matrixobj_bang(x);
}

static void matrix_size(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  int col, row;
  (void)s; /* unused */

  switch(argc) {
  case 0: /* size */
    if (x->m.row && x->m.col) {
      outlet_list(x->x_obj.ob_outlet, gensym("size"), 2, x->m.atombuffer);
    }
    break;
  case 1:
    row=atom_getfloat(argv);
    adjustsize(x, &x->m, row, row);
    matrix_set(&x->m, 0);
    break;
  default:
    row=atom_getfloat(argv++);
    col=atom_getfloat(argv);
    adjustsize(x, &x->m, row, col);
    matrix_set(&x->m, 0);
  }
}

void matrix_element(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  t_matrix*m = &x->m;
  t_atom *ap=m->atombuffer+2;
  int row=m->row, col=m->col;
  int r, c, i=row*col;
  (void)s; /* unused */

  switch (argc) {
  case 0:
    while(i--) {
      outlet_float(x->x_obj.ob_outlet, atom_getfloat(ap++));
    }
    break;
  case 1:
    r=c=atom_getfloat(argv)-1;
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    outlet_float(x->x_obj.ob_outlet, atom_getfloat(m->atombuffer+2+c+r*col));
    break;
  case 2:
    r=atom_getfloat(argv++)-1;
    c=atom_getfloat(argv++)-1;
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    outlet_float(x->x_obj.ob_outlet, atom_getfloat(m->atombuffer+2+c+r*col));
    break;
  default:
    r=atom_getfloat(argv++)-1;
    c=atom_getfloat(argv++)-1;
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    SETFLOAT(m->atombuffer+2+c+r*col, atom_getfloat(argv));
  }
}


void matrix_row(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  t_matrix*m = &x->m;
  t_atom *ap;
  int row=m->row, col=m->col;
  int r, c;
  t_float f;
  (void)s; /* unused */

  switch (argc) {
  case 0:
    /* row: get all rows as lists */
    for (r=0; r<row; r++) {
      outlet_list(x->x_obj.ob_outlet, gensym("row"), col,
                  m->atombuffer+r*col+2);
    }
    break;
  case 1:
    /* row <index>: get row<index> as list */
    r=atom_getfloat(argv)-1;
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    outlet_list(x->x_obj.ob_outlet, gensym("row"), col, m->atombuffer+r*col+2);
    break;
  case 2:
    /* row <index> <value>: set all elements of row<index> to <value> */
    r=atom_getfloat(argv)-1;
    f=atom_getfloat(argv+1);
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    for(c=0; c<col; c++) {
      SETFLOAT(m->atombuffer+2+c+r*col, f);
    }
    break;
  default:
    /* row <index> <value>...: set elements of row<index> to <value1> <value2> ... */
    r=atom_getfloat(argv++)-1;
    if (argc--<col) {
      pd_error(x, "matrix: sparse rows not yet supported : use [mtx_check]");
      return;
    }
    if ((r<0)||(r>=row)) {
      pd_error(x, "matrix: row index %d is out of range", r+1);
      return;
    }
    ap=m->atombuffer+2+col*r;
    memcpy(ap, argv, col*sizeof(t_atom));
    break;
  }
}

void matrix_col(t_matrixobj*x, t_symbol *s, int argc, t_atom *argv)
{
  t_matrix*m = &x->m;
  t_atom *ap;
  int row=m->row, col=m->col;
  int c, r;
  t_float f;
  (void)s; /* unused */

  switch (argc) {
  case 0:
    /* col: get all cols as lists */
    ap=(t_atom *)getbytes(row*sizeof(t_atom));
    for (c=0; c<col; c++) {
      for (r=0; r<row; r++) {
        SETFLOAT(&ap[r], atom_getfloat(m->atombuffer+2+c+col*r));
      }
      outlet_list(x->x_obj.ob_outlet, gensym("col"), row, ap);
    }
    freebytes(ap, row*sizeof(t_atom));
    break;
  case 1:
    /* col <index>: get col<index> as list */
    ap=(t_atom *)getbytes(row*sizeof(t_atom));
    c=atom_getfloat(argv)-1;
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    for (r=0; r<row; r++) {
      SETFLOAT(&ap[r], atom_getfloat(m->atombuffer+2+c+col*r));
    }
    outlet_list(x->x_obj.ob_outlet, gensym("col"), row, ap);
    freebytes(ap, row*sizeof(t_atom));
    break;
  case 2:
    /* row <index> <value>: set all elements of row<index> to <value> */
    c=atom_getint(argv)-1;
    f=atom_getfloat(argv+1);
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    for(r=0; r<row; r++) {
      SETFLOAT(m->atombuffer+2+c+r*col, f);
    }
    break;
  default:
    /* row <index> <value>...: set elements of row<index> to <value1> <value2> ... */
    c=atom_getfloat(argv++)-1;
    if (argc--<row) {
      pd_error(x, "matrix: sparse cols not yet supported : use [mtx_check]");
      return;
    }
    if ((c<0)||(c>=col)) {
      pd_error(x, "matrix: col index %d is out of range", c+1);
      return;
    }
    argv+=argc-1;
    if (argc>row) {
      argc=row;
    }
    while(argc--) {
      ap=m->atombuffer+2+c+col*argc;
      SETFLOAT(ap, atom_getfloat(argv--));
    }
  }
}

/* ------------- file I/O ------------------ */

static void matrix_read(t_matrixobj *x, t_symbol *filename)
{
  t_binbuf *bbuf = binbuf_new();
  t_atom *ap;
  int n;

  if (binbuf_read_via_path(bbuf, filename->s_name,
                           canvas_getdir(x->x_canvas)->s_name, 0)) {
    pd_error(x,"[matrix]: failed to read '%s'", filename->s_name);
  }

  ap=binbuf_getvec(bbuf);
  n =binbuf_getnatom(bbuf)-1;

  if ((ap->a_type == A_SYMBOL) &&
      (!strcmp(ap->a_w.w_symbol->s_name,"matrix")
       || !strcmp(ap->a_w.w_symbol->s_name,"#matrix")) ) {
    matrix_matrix2(x, &x->m, n, ap+1);
  }

  binbuf_free(bbuf);
}

static void matrix_write(t_matrixobj *x, t_symbol *filename)
{
  t_atom *ap=x->m.atombuffer+2;
  char filnam[MAXPDSTRING];
  int rows = x->m.row, cols = x->m.col;
  FILE *f=0;

  sys_bashfilename(filename->s_name, filnam);

  /* open file */
  if (!(f = fopen(filnam, "w"))) {
    pd_error(x, "[matrix]: failed to open '%s'", filnam);
  } else {
    char *text=(char *)getbytes(sizeof(char)*MAXPDSTRING);
    int textlen;

    /* header:
     * we now write "#matrix" instead of "matrix",
     * so that these files can easily read by other
     * applications such as octave
     */
    snprintf(text, MAXPDSTRING, "#matrix %d %d\n", rows, cols);
    text[MAXPDSTRING-1]=0;
    textlen = strlen(text);
    if (fwrite(text, textlen*sizeof(char), 1, f) < 1) {
      pd_error(x, "[matrix]: failed to write '%s'", filnam);
      goto end;
    }

    while(rows--) {
      int c = cols;
      while (c--) {
        t_float val = atom_getfloat(ap++);
        snprintf(text, MAXPDSTRING, "%.15f ", val);
        text[MAXPDSTRING-1]=0;
        textlen=strlen(text);
        if (fwrite(text, textlen*sizeof(char), 1, f) < 1) {
          pd_error(x, "[matrix]: failed to write '%s'", filnam);
          goto end;
        }
      }
      if (fwrite("\n", sizeof(char), 1, f) < 1) {
        pd_error(x, "[matrix]: failed to write '%s'", filnam);
        goto end;
      }
    }
    freebytes(text, sizeof(char)*MAXPDSTRING);
  }

end:
  /* close file */
  if (f) {
    fclose(f);
  }
}

static void matrix_list(t_matrixobj *x, t_symbol *s, int argc, t_atom *argv)
{
  /* like matrix, but without col/row information, so the previous size is kept */
  int row=x->m.row, col=x->m.col;
  (void)s; /* unused */

  if(!(row && col)) {
    pd_error(x, "[matrix]: unknown matrix dimensions");
    return;
  }
  if (argc<row*col) {
    pd_error(x, "[matrix]: sparse matrices not yet supported : use [mtx_check]!");
    return;
  }

  memcpy(x->m.atombuffer+2, argv, row*col*sizeof(t_atom));
  matrixobj_bang(x);
}

static void *matrix_new(t_symbol *s, int argc, t_atom *argv)
{
  t_matrixobj *x = (t_matrixobj *)pd_new(matrix_class);
  int row, col;
  (void)s; /* unused */

  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("matrix"), gensym(""));
  outlet_new(&x->x_obj, 0);

  x->m.atombuffer   = 0;
  x->x_canvas = canvas_getcurrent();

  switch (argc) {
  case 0:
    row = col = 0;
    break;
  case 1:
    if (argv->a_type == A_SYMBOL) {
      matrix_read(x, argv->a_w.w_symbol);
      return(x);
    }
    row = col = atom_getfloat(argv);
    break;
  default:
    row = atom_getfloat(argv++);
    col = atom_getfloat(argv++);
  }

  if(row && col) {
    adjustsize(x, &x->m, row, col);
    matrix_set(&x->m, 0);
  }

  return (x);
}

void matrix_setup(void)
{
  matrix_class = class_new(gensym("matrix"), (t_newmethod)matrix_new,
                           (t_method)matrixobj_free, sizeof(t_matrixobj), 0, A_GIMME, 0);
  class_addcreator((t_newmethod)matrix_new, gensym("mtx"), A_GIMME, 0);
  class_addcreator((t_newmethod)matrix_new, gensym("iemmatrix"), A_GIMME, 0);

  /* the core : functions for matrices */
  class_addmethod  (matrix_class, (t_method)matrix_matrix, gensym("matrix"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_matrix2, gensym(""),
                    A_GIMME, 0);

  /* the basics : functions for creation */
  class_addmethod  (matrix_class, (t_method)matrix_size, gensym("size"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_eye, gensym("eye"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_diag, gensym("diag"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_ones, gensym("ones"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_zeros, gensym("zeros"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_egg, gensym("egg"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrixobj_diegg, gensym("diegg"),
                    A_GIMME, 0);

  /* the rest : functions for anything */
  class_addbang    (matrix_class, matrixobj_bang);
  class_addfloat   (matrix_class, matrix_float);
  class_addlist    (matrix_class, matrix_list);
  class_addmethod  (matrix_class, (t_method)matrix_row, gensym("row"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrix_col, gensym("column"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrix_col, gensym("col"),
                    A_GIMME, 0);
  class_addmethod  (matrix_class, (t_method)matrix_element,
                    gensym("element"), A_GIMME, 0);

  /* the file functions */
  class_addmethod  (matrix_class, (t_method)matrix_write, gensym("write"),
                    A_SYMBOL, 0);
  class_addmethod  (matrix_class, (t_method)matrix_read , gensym("read") ,
                    A_SYMBOL, 0);
}

void iemtx_matrix_setup(void)
{
  matrix_setup();
}

void mtx_matrix_setup(void)
{
  matrix_setup();
}

void iematrix_setup(void)
{
  matrix_setup();
}
