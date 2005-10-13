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

static t_class *mtx_colon_class;
static t_symbol *col_sym;
static t_symbol *col_sym2;

typedef struct _MTXColon_ MTXColon;
struct _MTXColon_
{
   t_object x_obj;
   int size;
   t_symbol *colon_mode;

   t_atom *list_out;
   t_outlet *list_outlet;
};

static void deleteMTXColon (MTXColon *mtx_colon_obj) 
{
   if (mtx_colon_obj->list_out)
      freebytes (mtx_colon_obj->list_out, sizeof(t_atom)*(mtx_colon_obj->size+2));
}

static void mTXSetColonMode (MTXColon *mtx_colon_obj, t_symbol *c_mode) 
{
   mtx_colon_obj->colon_mode = c_mode;
}

static void *newMTXColon (t_symbol *s, int argc, t_atom *argv)
{
   MTXColon *mtx_colon_obj = (MTXColon *) pd_new (mtx_colon_class);

   mtx_colon_obj->list_outlet = outlet_new (&mtx_colon_obj->x_obj, gensym("matrix"));
   return ((void *) mtx_colon_obj);
} 

static void mTXColonBang (MTXColon *mtx_colon_obj)
{
   if (mtx_colon_obj->list_out) 
      outlet_anything(mtx_colon_obj->list_outlet, gensym("matrix"), 
	    mtx_colon_obj->size+2, mtx_colon_obj->list_out);
}

static void mTXColonList (MTXColon *mtx_colon_obj, t_symbol *s,
      int argc, t_atom *argv)
{
   int size;
   t_float startval;
   t_float stopval;
   t_float step;
   t_atom *list_out = mtx_colon_obj->list_out;
   if (argc == 3) {
      startval = atom_getfloat(argv++);
      step = atom_getfloat(argv++);
      stopval = atom_getfloat(argv++);
   }
   else if (argc == 2) {
      startval = atom_getfloat(argv++);
      stopval = atom_getfloat(argv++);
      step = 1.0f;
   }
   else {
      post("mtx_colon: wrong number of input arguments");
      return;
   }
      
   size = (int)((stopval- startval + step) / step);
   //post("startval %f stopval %f step %f, size %d",startval, stopval, step, size);
   if (size) {
      if (size!=mtx_colon_obj->size) {
	 if (list_out) 
	    list_out = (t_atom *) resizebytes (list_out,
		  sizeof(t_atom)*(mtx_colon_obj->size+2),
		  sizeof(t_atom)*(size+2));
	 else
	    list_out = (t_atom*) getbytes (sizeof(t_atom)*(size+2));
	 mtx_colon_obj->size = size;
      }
      mtx_colon_obj->list_out = list_out;

      if ((mtx_colon_obj->colon_mode == col_sym)||
	    (mtx_colon_obj->colon_mode == col_sym2)) {
	 SETFLOAT (&list_out[1],1.0f);
	 SETFLOAT (&list_out[0],(t_float)size);
      }
      else {
	 SETFLOAT (&list_out[0],1.0f);
	 SETFLOAT (&list_out[1],(t_float)size);
      }
      list_out += 2;
      for (;size--;list_out++,startval+=step)
	 SETFLOAT(list_out,startval);

      mTXColonBang (mtx_colon_obj);
   }
}

void mtx_colon_setup (void)
{
   mtx_colon_class = class_new 
      (gensym("mtx_colon"),
       (t_newmethod) newMTXColon,
       (t_method) deleteMTXColon,
       sizeof (MTXColon),
       CLASS_DEFAULT, A_GIMME, 0);
   class_addbang (mtx_colon_class, (t_method) mTXColonBang);
   class_addmethod (mtx_colon_class, (t_method) mTXSetColonMode, gensym("mode"), A_DEFSYMBOL, 0);
   class_addlist (mtx_colon_class, (t_method) mTXColonList);
   class_addcreator ((t_newmethod) newMTXColon, gensym("mtx_:"), A_GIMME, 0);
   class_sethelpsymbol (mtx_colon_class, gensym("iemmatrix/mtx_colon"));
   col_sym = gensym("col");
   col_sym2 = gensym("column");
}

void iemtx_colon_setup(void){
  mtx_colon_setup();
}