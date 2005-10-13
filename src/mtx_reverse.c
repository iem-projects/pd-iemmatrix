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

static t_class *mtx_reverse_class;
static t_symbol *row_sym;
static t_symbol *col_sym;
static t_symbol *col_sym2;

typedef struct _MTXreverse_ MTXreverse;
struct _MTXreverse_
{
   t_object x_obj;
   int size;
   //int reverse_dimension;
   t_symbol *reverse_mode;

   t_outlet *list_outlet;

   t_atom *list_out;
   t_atom *list_in;
};

static void deleteMTXreverse (MTXreverse *mtx_reverse_obj) 
{
   if (mtx_reverse_obj->list_out)
      freebytes (mtx_reverse_obj->list_out, sizeof(t_atom)*(mtx_reverse_obj->size+2));
}
static void mTXSetReverseMode (MTXreverse *mtx_reverse_obj, t_symbol *c_mode)
{
   mtx_reverse_obj->reverse_mode = c_mode;
}
/*
static void mTXSetreverseDimension (MTXreverse *mtx_reverse_obj, t_float c_dim)
{
   int dimension = (int) c_dim;
   dimension = (dimension > 0)?dimension:0;
   dimension = (dimension < 2)?dimension:2;
   mtx_reverse_obj->reverse_dimension = dimension;
}
*/


static void *newMTXreverse (t_symbol *s, int argc, t_atom *argv)
{
   MTXreverse *mtx_reverse_obj = (MTXreverse *) pd_new (mtx_reverse_class);
   mTXSetReverseMode (mtx_reverse_obj, gensym(":"));
   switch ((argc>1)?1:argc) {
      case 1:
	 mTXSetReverseMode (mtx_reverse_obj, atom_getsymbol (argv));
   }
   /*int c_dim = 0;

   mtx_reverse_obj->reverse_dimension = c_dim;
   switch ((argc>1)?1:argc) {
      case 1:
	 c_dim = atom_getint(argv);
   }
   mTXSetreverseDimension (mtx_reverse_obj, (t_float) c_dim);
*/
   mtx_reverse_obj->list_outlet = outlet_new (&mtx_reverse_obj->x_obj, gensym("matrix"));
   return ((void *) mtx_reverse_obj);
} 

static void mTXreverseBang (MTXreverse *mtx_reverse_obj)
{
   if (mtx_reverse_obj->list_out) 
      outlet_anything(mtx_reverse_obj->list_outlet, gensym("matrix"), 
	    mtx_reverse_obj->size+2, mtx_reverse_obj->list_out);
}

static void copyList (int n, t_atom *x, t_atom *y)
{
   for (;n--;)  
      *y++ = *x++;
}
static void reverseList (int n, t_atom *y)
{
   t_atom *read = y;
   t_atom tmp;
   y += n-1;
   n >>= 1;
   for (;n--;) { 
      tmp = *y;
      *y-- = *read;
      *read++ = tmp;
   }
}
static void reverseListStep (int n, int step, t_atom *y)
{
   t_atom *read = y;
   t_atom tmp;
   n /= step;
   y += (n-1) * step;
   n >>= 1;
   for (;n--; y-=step, read+=step) { 
      tmp = *y;
      *y = *read;
      *read = tmp;
   }
}

static void mTXreverseMatrix (MTXreverse *mtx_reverse_obj, t_symbol *s, 
      int argc, t_atom *argv)
{
   int rows = atom_getint (argv++);
   int columns = atom_getint (argv++);
   int size = rows * columns;
   int list_size = argc - 2;
   t_atom *list_in = argv;
   t_atom *list_out = mtx_reverse_obj->list_out;
   int count;

   // size check
   if (!size) {
      post("mtx_reverse: invalid dimensions");
      return;
   }
   else if (list_size<size) {
      post("mtx_reverse: sparse matrix not yet supported: use \"mtx_check\"");
      return;
   }
   
   if (size != mtx_reverse_obj->size) {
      if (!list_out)
	 list_out = (t_atom *) getbytes (sizeof (t_atom) * (size + 2));
      else
	 list_out = (t_atom *) resizebytes (list_out,
	       sizeof (t_atom) * (mtx_reverse_obj->size+2),
	       sizeof (t_atom) * (size + 2));
   }

   mtx_reverse_obj->size = size;
   mtx_reverse_obj->list_out = list_out;

   // main part
   list_out += 2;
   copyList (size, argv, list_out);

   if ((mtx_reverse_obj->reverse_mode == col_sym)||
	 (mtx_reverse_obj->reverse_mode == col_sym2)) {
      for (count = columns; count--; list_out++)
	 reverseListStep (size, columns, list_out);
   }
   else if (mtx_reverse_obj->reverse_mode == row_sym) {
      for (count = rows; count--; list_out += columns) 
	 reverseList (columns, list_out);
   }
   else 
      reverseList (size, list_out); 

/*
   switch (mtx_reverse_obj->reverse_dimension) {
      case 2:
	 for (count = columns; count--; list_out++)
	    reverseListStep (size, columns, list_out);
	 break;
      case 1:
	 for (count = rows; count--; list_out += columns) 
	    reverseList (columns, list_out);
	 break;
      case 0:
	 reverseList (size, list_out); 
	 break;
   }
   */
   list_out = mtx_reverse_obj->list_out;


   SETSYMBOL(list_out, gensym("matrix"));
   SETFLOAT(list_out, rows);
   SETFLOAT(&list_out[1], columns);
   outlet_anything(mtx_reverse_obj->list_outlet, gensym("matrix"), 
	 mtx_reverse_obj->size+2, list_out);
}

void mtx_reverse_setup (void)
{
   mtx_reverse_class = class_new 
      (gensym("mtx_reverse"),
       (t_newmethod) newMTXreverse,
       (t_method) deleteMTXreverse,
       sizeof (MTXreverse),
       CLASS_DEFAULT, A_GIMME, 0);
   class_addbang (mtx_reverse_class, (t_method) mTXreverseBang);
   class_addmethod (mtx_reverse_class, (t_method) mTXreverseMatrix, gensym("matrix"), A_GIMME,0);
//   class_addmethod (mtx_reverse_class, (t_method) mTXSetreverseDimension, gensym("dimension"), A_DEFFLOAT,0);
   class_addmethod (mtx_reverse_class, (t_method) mTXSetReverseMode, gensym("mode"), A_DEFSYMBOL,0);
   class_sethelpsymbol (mtx_reverse_class, gensym("iemmatrix/mtx_reverse"));
   row_sym = gensym("row");
   col_sym = gensym("col");
   col_sym2 = gensym("column");
}

void iemtx_reverse_setup(void){
  mtx_reverse_setup();
}