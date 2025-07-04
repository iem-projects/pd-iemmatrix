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

static t_class *mtx_cumsum_class;
static t_symbol *row_sym;
static t_symbol *col_sym;
static t_symbol *col_sym2;

typedef struct _MTXCumsum_ MTXCumsum;
struct _MTXCumsum_ {
  t_object x_obj;
  int rows;
  int columns;
  int size;
  int cumsum_direction;
  t_symbol *cumsum_mode;

  t_outlet *list_outlet;

  t_atom *list_out;
  t_atom *list_in;
  t_float *x;
  t_float *y;
};

static void deleteMTXCumsum (MTXCumsum *mtx_cumsum_obj)
{
  if (mtx_cumsum_obj->list_out) {
    freebytes (mtx_cumsum_obj->list_out,
               sizeof(t_atom)*(mtx_cumsum_obj->size+2));
  }
  if (mtx_cumsum_obj->x) {
    freebytes (mtx_cumsum_obj->x, sizeof(t_float)*(mtx_cumsum_obj->size));
  }
  if (mtx_cumsum_obj->y) {
    freebytes (mtx_cumsum_obj->y, sizeof(t_float)*(mtx_cumsum_obj->size));
  }
}

static void mTXSetCumsumDirection (MTXCumsum *mtx_cumsum_obj,
                                   t_float c_dir)
{
  int direction = (int) c_dir;
  mtx_cumsum_obj->cumsum_direction = (direction==-1)?direction:1;
}

static void mTXSetCumsumMode (MTXCumsum *mtx_cumsum_obj, t_symbol *m_sym)
{
  mtx_cumsum_obj->cumsum_mode = m_sym;
}

static void *newMTXCumsum (t_symbol *s, int argc, t_atom *argv)
{
  MTXCumsum *mtx_cumsum_obj = (MTXCumsum *) pd_new (mtx_cumsum_class);
  mTXSetCumsumMode (mtx_cumsum_obj, gensym(":"));
  mTXSetCumsumDirection (mtx_cumsum_obj, 1.0f);
  if (argc>=1) {
    if (argv[0].a_type == A_SYMBOL) {
      mTXSetCumsumMode (mtx_cumsum_obj, atom_getsymbol (argv));
      if (argc>=2) {
        if (argv[1].a_type != A_SYMBOL) {
          mTXSetCumsumDirection (mtx_cumsum_obj, atom_getfloat (argv+1));
        } else {
          pd_error(mtx_cumsum_obj, "[%s]: 2nd arg ignored. supposed to be float", s->s_name);
        }
      }
    } else {
      mTXSetCumsumDirection (mtx_cumsum_obj, atom_getfloat (argv));
      if (argc>=2) {
        if (argv[1].a_type == A_SYMBOL) {
          mTXSetCumsumMode (mtx_cumsum_obj, atom_getsymbol (argv+1));
        } else {
          pd_error(mtx_cumsum_obj, "[%s]: 2nd arg ignored. supposed to be symbolic, e.g. \"row\", \"col\", \":\"",
                   s->s_name);
        }
      }
    }
  }

  mtx_cumsum_obj->list_outlet = outlet_new (&mtx_cumsum_obj->x_obj,
                                gensym("matrix"));
  return ((void *) mtx_cumsum_obj);
}

static void mTXCumsumBang (MTXCumsum *mtx_cumsum_obj)
{
  if (mtx_cumsum_obj->list_out)
    outlet_anything(mtx_cumsum_obj->list_outlet, gensym("matrix"),
                    mtx_cumsum_obj->size+2, mtx_cumsum_obj->list_out);
}

static void cumSum (int n, t_float *x, t_float *y)
{
  t_float accu = 0.0f;
  for (; n--; x++, y++) {
    accu += *x;
    *y = accu;
  }
}
static void cumSumReverse (int n, t_float *x, t_float *y)
{
  t_float accu = 0.0f;
  for (; n--; x--, y--) {
    accu += *x;
    *y = accu;
  }
}

static void mTXCumsumMatrix (MTXCumsum *mtx_cumsum_obj, t_symbol *s,
                             int argc, t_atom *argv)
{
  int rows, columns, size;
  t_atom *list_ptr = argv+2;
  t_atom *list_out = mtx_cumsum_obj->list_out;
  t_float *x = mtx_cumsum_obj->x;
  t_float *y = mtx_cumsum_obj->y;
  int count;
  (void)s; /* unused */

  /* size check */
  if(iemmatrix_check(mtx_cumsum_obj, s, argc, argv, 0))return;
  rows = atom_getint (argv+0);
  columns = atom_getint (argv+1);
  size = rows * columns;

  if ((!x)||(!list_out)||(!y)) {
    if (!x) {
      x = (t_float *) getbytes (sizeof (t_float) * (size));
    }
    if (!y) {
      y = (t_float *) getbytes (sizeof (t_float) * (size));
    }
    if (!list_out) {
      list_out = (t_atom *) getbytes (sizeof (t_atom) * (size+2));
    }
  } else if (size != mtx_cumsum_obj->size) {
    x = (t_float *) resizebytes (x,
                                 sizeof (t_float) * (mtx_cumsum_obj->size),
                                 sizeof (t_float) * (size));
    y = (t_float *) resizebytes (y,
                                 sizeof (t_float) * (mtx_cumsum_obj->size),
                                 sizeof (t_float) * (size));
    list_out = (t_atom *) resizebytes (list_out,
                                       sizeof (t_atom) * (mtx_cumsum_obj->size+2),
                                       sizeof (t_atom) * (size + 2));
  }
  mtx_cumsum_obj->size = size;
  mtx_cumsum_obj->rows = rows;
  mtx_cumsum_obj->columns = columns;
  mtx_cumsum_obj->list_out = list_out;
  mtx_cumsum_obj->x = x;
  mtx_cumsum_obj->y = y;

  /* main part */
  /* reading matrix from inlet */
  if ((mtx_cumsum_obj->cumsum_mode == col_sym) ||
      (mtx_cumsum_obj->cumsum_mode == col_sym2)) {
    iemmatrix_list2floats_modulo(x, list_ptr, size, columns);
    columns = mtx_cumsum_obj->rows;
    rows = mtx_cumsum_obj->columns;
  } else {
    iemmatrix_list2floats(x, list_ptr, size);
  }

  /* calculating cumsum */
  if (mtx_cumsum_obj->cumsum_direction == -1) {
    if ((mtx_cumsum_obj->cumsum_mode == row_sym) ||
        (mtx_cumsum_obj->cumsum_mode == col_sym) ||
        (mtx_cumsum_obj->cumsum_mode == col_sym2)) {
      x += columns-1;
      y += columns-1;

      for (count = rows; count--; x += columns, y += columns) {
        cumSumReverse (columns,x,y);
      }
    } else {
      x += size-1;
      y += size-1;
      cumSumReverse (size, x, y);
    }
  } else if ((mtx_cumsum_obj->cumsum_mode == row_sym) ||
             (mtx_cumsum_obj->cumsum_mode == col_sym) ||
             (mtx_cumsum_obj->cumsum_mode == col_sym2))
    for (count = rows; count--; x += columns, y += columns) {
      cumSum (columns,x,y);
    }
  else {
    cumSum (size, x, y);
  }

  x = mtx_cumsum_obj->x;
  y = mtx_cumsum_obj->y;

  /* writing matrix to outlet */
  if ((mtx_cumsum_obj->cumsum_mode == col_sym) ||
      (mtx_cumsum_obj->cumsum_mode == col_sym2)) {
    columns = mtx_cumsum_obj->columns;
    rows = mtx_cumsum_obj->rows;
    iemmatrix_floats2list_modulo(list_out+2, y, size, columns);
  } else {
    iemmatrix_floats2list(list_out+2, y, size);
  }

  SETSYMBOL(list_out, gensym("matrix"));
  SETFLOAT(list_out, rows);
  SETFLOAT(&list_out[1], columns);
  outlet_anything(mtx_cumsum_obj->list_outlet, gensym("matrix"),
                  mtx_cumsum_obj->size+2, list_out);
}

void mtx_cumsum_setup (void)
{
  mtx_cumsum_class = class_new
                     (gensym("mtx_cumsum"),
                      (t_newmethod) newMTXCumsum,
                      (t_method) deleteMTXCumsum,
                      sizeof (MTXCumsum),
                      CLASS_DEFAULT, A_GIMME, 0);
  class_addbang (mtx_cumsum_class, (t_method) mTXCumsumBang);
  class_addmethod (mtx_cumsum_class, (t_method) mTXCumsumMatrix,
                   gensym("matrix"), A_GIMME,0);
  class_addmethod (mtx_cumsum_class, (t_method) mTXSetCumsumMode,
                   gensym("mode"), A_DEFSYMBOL,0);
  class_addmethod (mtx_cumsum_class, (t_method) mTXSetCumsumDirection,
                   gensym("direction"), A_DEFFLOAT,0);

  row_sym = gensym("row");
  col_sym = gensym("col");
  col_sym2 = gensym("column");
}

void iemtx_cumsum_setup(void)
{
  mtx_cumsum_setup();
}
