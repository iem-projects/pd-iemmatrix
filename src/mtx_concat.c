/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly referring to matlab/octave matrix functions
 *
 * Copyright (c) 2005, Franz Zotter
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#include "iemmatrix.h"

static t_class *mtx_concat_class;
static t_symbol *row_sym;

typedef struct _MTXconcat_ MTXconcat;
struct _MTXconcat_ {
  t_object x_obj;
  int size;
  int concat_mode;
  t_matrix mtx_in1;
  t_matrix mtx_in2;
  t_matrix mtx_out;

  t_outlet *outl;
};

static void deleteMTXConcat (MTXconcat *mtx_concat_obj)
{
  matrix_free(&mtx_concat_obj->mtx_in2);
  matrix_free(&mtx_concat_obj->mtx_out);
}

static void mTXSetConcatMode (MTXconcat *mtx_concat_obj, t_symbol *c_mode)
{
  char c=*c_mode->s_name;
  switch(c) {
  case 'c':
  case 'C':
  case ':': /* "column" */
    mtx_concat_obj->concat_mode = 1;
    break;
  case 'r':
  case 'R': /* "row" */
    mtx_concat_obj->concat_mode = 0;
    break;
  default:
    pd_error(mtx_concat_obj, "mtx_concat: invalid mode '%s'", c_mode->s_name);
    break;
  }
}

static void *newMTXConcat (t_symbol *s, int argc, t_atom *argv)
{
  MTXconcat *mtx_concat_obj = (MTXconcat *) pd_new (mtx_concat_class);
  (void)s; /* unused */
  if(argc&&(A_SYMBOL==argv->a_type)) {
    mTXSetConcatMode (mtx_concat_obj, atom_getsymbol (argv));
  } else {
    mTXSetConcatMode (mtx_concat_obj, gensym(":"));
  }

  mtx_concat_obj->outl = mtx_concat_obj->mtx_out.x_outlet = outlet_new (
                           &mtx_concat_obj->x_obj, gensym("matrix"));
  inlet_new(&mtx_concat_obj->x_obj, &mtx_concat_obj->x_obj.ob_pd,
            gensym("matrix"),gensym(""));
  return ((void *) mtx_concat_obj);
}

static void mTXConcatBang (MTXconcat *mtx_concat_obj)
{
  outlet_anything(mtx_concat_obj->outl, gensym("matrix"),
                  mtx_concat_obj->mtx_out.row * mtx_concat_obj->mtx_out.col + 2,
                  mtx_concat_obj->mtx_out.atombuffer);
}

static void mTXConcatMatrix2 (MTXconcat *mtx_concat_obj, t_symbol *s,
                              int argc, t_atom *argv)
{
  matrix_matrix2 (&mtx_concat_obj->mtx_in2, s, argc, argv);
}

static void mTXConcatDoRowConcatenation (MTXconcat *mtx_concat_obj,
    t_matrix *mtx1, t_matrix *mtx2, t_matrix *mtx_out)
{
  int mcols = mtx1->col + mtx2->col;
  int cnt;
  t_atom *ptr_in1 = mtx1->atombuffer+2;
  t_atom *ptr_in2 = mtx2->atombuffer+2;
  t_atom *ptr_out;

  if (mtx1->row != mtx2->row) {
    pd_error(mtx_concat_obj, "[mtx_concat]: row-mode: matrices must have same number of rows!");
    return;
  }
  adjustsize (mtx_out, mtx1->row, mcols);
  ptr_out = mtx_out->atombuffer+2;
  for (cnt=mtx1->row; cnt--; ptr_in1 += mtx1->col,
       ptr_in2 += mtx2->col, ptr_out += mtx_out->col) {
    memcpy (ptr_out, ptr_in1, mtx1->col * sizeof(t_atom));
    memcpy (ptr_out+mtx1->col, ptr_in2, mtx2->col * sizeof(t_atom));
  }
  mTXConcatBang(mtx_concat_obj);
}
static void mTXConcatDoColConcatenation (MTXconcat *mtx_concat_obj,
    t_matrix *mtx1, t_matrix *mtx2, t_matrix *mtx_out)
{
  int mrows = mtx1->row + mtx2->row;
  int cnt;
  t_atom *ptr_in1 = mtx1->atombuffer+2;
  t_atom *ptr_in2 = mtx2->atombuffer+2;
  t_atom *ptr_out;

  if (mtx1->col != mtx2->col) {
    pd_error(mtx_concat_obj, "[mtx_concat]: col-mode: matrices must have same number of columns!");
    return;
  }
  adjustsize (mtx_out, mrows, mtx1->col);
  ptr_out = mtx_out->atombuffer+2;
  for (cnt=mtx1->row; cnt--; ptr_in1 += mtx1->col,
       ptr_out += mtx_out->col) {
    memcpy (ptr_out, ptr_in1, mtx1->col * sizeof(t_atom));
  }
  for (cnt=mtx2->row; cnt--; ptr_in2 += mtx2->col,
       ptr_out += mtx_out->col) {
    memcpy (ptr_out, ptr_in2, mtx2->col * sizeof(t_atom));
  }

  mTXConcatBang(mtx_concat_obj);
}
static void mTXConcatMatrix (MTXconcat *mtx_concat_obj, t_symbol *s,
                             int argc, t_atom *argv)
{
  int rows = atom_getint (argv);
  int columns = atom_getint (argv+1);
  t_matrix *mtx_in1 = &mtx_concat_obj->mtx_in1;
  t_matrix *mtx_in2 = &mtx_concat_obj->mtx_in2;
  t_matrix *mtx_out = &mtx_concat_obj->mtx_out;
  (void)s; /* unused */

  /* size check */
  if(iemmatrix_check(mtx_concat_obj, s, argc, argv, 0))return;

  mtx_in1->row = rows;
  mtx_in1->col = columns;
  mtx_in1->atombuffer = argv;
  /* alternatively to the above: */
  /*   matrix_matrix2 (mtx_in1, s, argc, argv); */

  if (mtx_concat_obj->concat_mode == 0) {
    mTXConcatDoRowConcatenation(mtx_concat_obj, mtx_in1, mtx_in2, mtx_out);
  } else {
    mTXConcatDoColConcatenation(mtx_concat_obj, mtx_in1, mtx_in2, mtx_out);
  }
}

void mtx_concat_setup (void)
{
  mtx_concat_class = class_new
                     (gensym("mtx_concat"),
                      (t_newmethod) newMTXConcat,
                      (t_method) deleteMTXConcat,
                      sizeof (MTXconcat),
                      CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_concat_class, (t_method) mTXConcatBang);
  class_addmethod (mtx_concat_class, (t_method) mTXConcatMatrix,
                   gensym("matrix"), A_GIMME,0);
  class_addmethod (mtx_concat_class, (t_method) mTXConcatMatrix2, gensym(""),
                   A_GIMME,0);
  class_addmethod (mtx_concat_class, (t_method) mTXSetConcatMode,
                   gensym("mode"), A_DEFSYMBOL,0);

  row_sym = gensym("row");
}

void iemtx_concat_setup(void)
{
  mtx_concat_setup();
}
