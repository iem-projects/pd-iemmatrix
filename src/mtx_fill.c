/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  mostly refering to matlab/octave matrix functions
 *
 * Copyright (c) 2005, Franz Zotter
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */


#include "iemmatrix.h"

typedef enum {
   FILL_SUBMATRIX,
   FILL_INDEXED_ELEMENTS
} FillStyle;

static t_class *mtx_fill_class;

typedef struct _MTXfill_ MTXfill;
struct _MTXfill_
{
   t_object x_obj;
   int size;
   int rows;
   int columns;
   
   int fill_startcol;
   int fill_startrow;

   int *index;
   int index_size;
   int max_index;

   FillStyle fill_type;

   t_outlet *list_outlet;

   t_atom *list_in;
   t_atom *list_out;
};

static void deleteMTXFill (MTXfill *mtx_fill_obj) 
{
   if (mtx_fill_obj->list_in)
      freebytes (mtx_fill_obj->list_in, sizeof(t_atom)*(mtx_fill_obj->size+2));
   if (mtx_fill_obj->list_out)
      freebytes (mtx_fill_obj->list_out, sizeof(t_atom)*(mtx_fill_obj->size+2));
   if (mtx_fill_obj->index)
      freebytes (mtx_fill_obj->index, sizeof(int)*(mtx_fill_obj->index_size));
}

static void setListConstFloat (int size, t_float f, t_atom *y)
{
   for(;size--;y++)
      SETFLOAT(y,f);
}
static void copyList (int size, t_atom *x, t_atom *y)
{
   while(size--)
      *y++=*x++;
}
static int copyNonZeroAtomsToIntegerArrayMax (int *size, t_atom *x, int *y)
{
   int idx;
   int n = *size;
   int max = atom_getint(x);
   *size = 0;
   for (;n--;x++) {
      idx = atom_getint (x);
      if (idx) {
	 size[0]++;
	 *y++ = idx;
         max = (idx > max)?idx:max;
      }
   }
   return max;
}

static void writeIndexedValuesIntoMatrix (int n, int *index, t_atom *x, t_atom *y)
{
   for (;n--;index++,x++)
      if (*index)
	 y[*index-1] = *x;
}
static void writeFloatIndexedIntoMatrix (int n, int *index, t_float f, t_atom *y)
{
   for (;n--;index++)
      if (*index)
	 SETFLOAT(&y[*index-1], f);
}

static void mTXFillIndexMatrix (MTXfill *mtx_fill_obj, t_symbol *s, 
      int argc, t_atom *argv)
{
   int rows = atom_getint (argv++);
   int columns = atom_getint (argv++);
   int size = rows * columns;
   int list_size = argc - 2;
   int *index = mtx_fill_obj->index;

   // size check
   if (!size) {
      post("mtx_fill: invalid dimensions/invalid start index");
      return;
   }
   
   if (list_size == 0) {
      if ((rows<1) || (columns<1)){
	 post("mtx_fill: row and column indices must be >0");
	 return;
      }
      mtx_fill_obj->fill_startrow = rows;
      mtx_fill_obj->fill_startcol = columns;
      mtx_fill_obj->fill_type = FILL_SUBMATRIX;
   }
   else if (list_size<size) {
      post("mtx_fill: sparse matrix not yet supported: use \"mtx_check\"");
      return;
   }
   else {
      if (size > mtx_fill_obj->index_size) {
	 if (!index)
	    index = (int *) getbytes (sizeof (int) * (size + 2));
	 else
	    index = (int *) resizebytes (index,
		  sizeof (int) * (mtx_fill_obj->index_size+2),
		  sizeof (t_atom) * (size + 2));
	 mtx_fill_obj->index_size = size;
      }
      mtx_fill_obj->max_index = 
	 copyNonZeroAtomsToIntegerArrayMax (&size, argv++, index);
      if (!size) {
	 post("mtx_fill: indexing matrix contains zero-values only!!!");
	 return;
      }
      if (size != mtx_fill_obj->index_size) {
	 index = (int *)  resizebytes (index,
		  sizeof (int) * (mtx_fill_obj->index_size+2),
		  sizeof (t_atom) * (size + 2));
	 mtx_fill_obj->index_size = size;
      }
      mtx_fill_obj->fill_type = FILL_INDEXED_ELEMENTS;
      mtx_fill_obj->index = index;
   }
}

static void *newMTXFill (t_symbol *s, int argc, t_atom *argv)
{
   MTXfill *mtx_fill_obj = (MTXfill *) pd_new (mtx_fill_class);
  
   mtx_fill_obj->fill_startrow = 1;
   mtx_fill_obj->fill_startcol = 1;
   mtx_fill_obj->fill_type = FILL_SUBMATRIX;
   error("[mtx_fill]: this object _might_ change in the future!");
   if (argc) {
      if (atom_getsymbol(argv)==gensym("matrix")) 
	 mTXFillIndexMatrix (mtx_fill_obj, s, argc-1, argv+1);
      else
	 pd_error(mtx_fill_obj, "mtx_fill: creation argument must be 'matrix <startrow> <startcol>' for submatrix filling or 'matrix rows columns [...]' for indexed filling with scalar/matrices"); 
   }

   mtx_fill_obj->list_outlet = outlet_new (&mtx_fill_obj->x_obj, gensym("matrix"));
   inlet_new(&mtx_fill_obj->x_obj, &mtx_fill_obj->x_obj.ob_pd, gensym("matrix"),gensym("fill_mtx"));
   inlet_new(&mtx_fill_obj->x_obj, &mtx_fill_obj->x_obj.ob_pd, gensym("matrix"),gensym("index"));
   return ((void *) mtx_fill_obj);
} 

static void mTXBigMatrix (MTXfill *mtx_fill_obj, t_symbol *s, 
      int argc, t_atom *argv)
{
   int rows = atom_getint (argv++);
   int columns = atom_getint (argv++);
   int size = rows * columns;
   int list_size = argc - 2;
   t_atom *list_in = mtx_fill_obj->list_in;
   t_atom *list_out = mtx_fill_obj->list_out;

   // size check
   if (!size) {
      post("mtx_fill: invalid dimensions");
      return;
   }
   else if (list_size<size) {
      post("mtx_fill: sparse matrix not yet supported: use \"mtx_check\"");
      return;
   }
   
   if (size != mtx_fill_obj->size) {
      if (!list_out)
	 list_out = (t_atom *) getbytes (sizeof (t_atom) * (size + 2));
      else
	 list_out = (t_atom *) resizebytes (list_out,
	       sizeof (t_atom) * (mtx_fill_obj->size+2),
	       sizeof (t_atom) * (size + 2));
      if (!list_in)
	 list_in = (t_atom *) getbytes (sizeof (t_atom) * (size + 2));
      else
	 list_in = (t_atom *) resizebytes (list_in,
	       sizeof (t_atom) * (mtx_fill_obj->size+2),
	       sizeof (t_atom) * (size + 2));
   }

   mtx_fill_obj->size = size;
   mtx_fill_obj->columns = columns;
   mtx_fill_obj->rows = rows;
   mtx_fill_obj->list_out = list_out;
   mtx_fill_obj->list_in = list_in;

   copyList (size, argv, list_in);
}


static void mTXFillBang (MTXfill *mtx_fill_obj)
{
   if (mtx_fill_obj->list_out) 
      outlet_anything(mtx_fill_obj->list_outlet, gensym("matrix"), 
	    mtx_fill_obj->size+2, mtx_fill_obj->list_out);
}

static void writeFillMatrixIntoList (int fillrows, const int fillcols, int columns, t_atom *x, t_atom *y)
{
   for (;fillrows--;x+=fillcols,y+=columns)
      copyList(fillcols, x, y);
}

static void mTXFillScalar (MTXfill *mtx_fill_obj, t_float f)
{
   t_atom *list_out = mtx_fill_obj->list_out;
   t_atom *list_in = mtx_fill_obj->list_in;
   int rows = mtx_fill_obj->rows;
   int columns = mtx_fill_obj->columns;
   if (mtx_fill_obj->fill_type == FILL_INDEXED_ELEMENTS) {
      if (mtx_fill_obj->max_index > mtx_fill_obj->size) {
	 post("mtx_fill: index matrix index exceeds matrix borders");
	 return;
      }
      else if (mtx_fill_obj->size == 0) {
	 post("mtx_fill: no matrix defined for filling");
	 return;
      }

      // main part
      list_out += 2;
      copyList (mtx_fill_obj->size, list_in, list_out);

      writeFloatIndexedIntoMatrix (mtx_fill_obj->index_size,
	    mtx_fill_obj->index, f,list_out);
      list_out = mtx_fill_obj->list_out;
      SETSYMBOL(list_out, gensym("matrix"));
      SETFLOAT(list_out, rows);
      SETFLOAT(&list_out[1], columns);
      outlet_anything(mtx_fill_obj->list_outlet, gensym("matrix"), 
	    mtx_fill_obj->size+2, list_out);
   }
   else
      post("mtx_fill: scalar fill for submatrices not supported yet");
}


static void mTXFillMatrix (MTXfill *mtx_fill_obj, t_symbol *s, 
      int argc, t_atom *argv)
{
   int fill_rows = atom_getint (argv++);
   int fill_columns = atom_getint (argv++);
   int fill_size = fill_rows * fill_columns;
   int list_size = argc - 2;
   int rows = mtx_fill_obj->rows;
   int columns = mtx_fill_obj->columns;
   t_atom *fill_mtx = argv;
   t_atom *list_in = mtx_fill_obj->list_in;
   t_atom *list_out = mtx_fill_obj->list_out;
   int stopcol = mtx_fill_obj->fill_startcol+fill_columns-1;
   int stoprow = mtx_fill_obj->fill_startrow+fill_rows-1;

   // size check
   if (!list_size) {
      post("mtx_fill: invalid dimensions");
      return;
   }
   switch (mtx_fill_obj->fill_type) {
      case FILL_SUBMATRIX:
	 if (list_size < fill_size) {
	    post("mtx_fill: sparse matrix not yet supported: use \"mtx_check\"");
	    return;
	 }
	 if ((stopcol > columns) ||
	       (stoprow > rows)) {
	    post("mtx_fill: fill matrix index exceeds matrix borders");
	    return;
	 }
	 break;
      case FILL_INDEXED_ELEMENTS:
	 if (list_size > mtx_fill_obj->index_size) {
	    post("mtx_fill: fill matrix smaller than indexing vector");
	    return;
	 }
	 else if (mtx_fill_obj->max_index > mtx_fill_obj->size) {
	    post("mtx_fill: index matrix index exceeds matrix borders");
	    return;
	 }
	 break;
   }
   if (mtx_fill_obj->size == 0) {
      post("mtx_fill: no matrix defined for filling");
      return;
   }
   

   // main part
   list_out += 2;
   copyList (mtx_fill_obj->size, list_in, list_out);

   switch (mtx_fill_obj->fill_type) {
      case FILL_SUBMATRIX:
	 list_out += columns * (mtx_fill_obj->fill_startrow-1) + 
	    mtx_fill_obj->fill_startcol-1;
	 writeFillMatrixIntoList (fill_rows, fill_columns, 
	       columns, fill_mtx, list_out);
	 break;
      case FILL_INDEXED_ELEMENTS:
	 writeIndexedValuesIntoMatrix (mtx_fill_obj->index_size,
	       mtx_fill_obj->index, fill_mtx,list_out);
	 break;
   }
   list_out = mtx_fill_obj->list_out;
   SETSYMBOL(list_out, gensym("matrix"));
   SETFLOAT(list_out, rows);
   SETFLOAT(&list_out[1], columns);
   outlet_anything(mtx_fill_obj->list_outlet, gensym("matrix"), 
	 mtx_fill_obj->size+2, list_out);
}

void mtx_fill_setup (void)
{
   mtx_fill_class = class_new 
      (gensym("mtx_fill"),
       (t_newmethod) newMTXFill,
       (t_method) deleteMTXFill,
       sizeof (MTXfill),
       CLASS_DEFAULT, A_GIMME, 0);
   class_addbang (mtx_fill_class, (t_method) mTXFillBang);
   class_addmethod (mtx_fill_class, (t_method) mTXFillMatrix, gensym("matrix"), A_GIMME,0);
   class_addmethod (mtx_fill_class, (t_method) mTXBigMatrix, gensym("fill_mtx"), A_GIMME,0);
   class_addmethod (mtx_fill_class, (t_method) mTXFillIndexMatrix, gensym("index"), A_GIMME,0);
   class_addfloat (mtx_fill_class, (t_method) mTXFillScalar);
   class_sethelpsymbol (mtx_fill_class, gensym("iemmatrix/mtx_fill"));
}

void iemtx_fill_setup(void){
  mtx_fill_setup();
}
