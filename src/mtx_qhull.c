/*
 *  iemmatrix
 *
 *  objects for manipulating simple matrices
 *  own qhull algorithm implementation
 *
 * Copyright (c) 2012, Franz Zotter
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */


#include "iemmatrix.h"
#include <stdlib.h>
#include "mtx_qhull/list.h"
#include "mtx_qhull/vectors.h"
#include "mtx_qhull/zhull.h"
#include "mtx_qhull/list.c"
#include "mtx_qhull/vectors.c"
#include "mtx_qhull/zhull.c"

static t_class *mtx_qhull_class;

typedef struct MTXQhull_ MTXQhull;
struct MTXQhull_  {
    t_object x_obj;
    t_outlet *outl;
    t_atom *list;
    size_t size;
    size_t hull_size;
    int iter;
    zhull_t *zh;
};

static void deleteMTXQhull(MTXQhull *xo) {
    if (xo->zh!=0) {
        free(xo->zh);
        xo->zh=0;
    }
    if (xo->list!=0) {
        free(xo->list);
        xo->list=0;
        xo->size=0;
    }
}

void mTXQhullSetIterations(MTXQhull *xo, t_float f) {
    xo->iter=(int)f;
    if (xo->iter<0)
        xo->iter=0;
}

static void *newMTXQhull(t_symbol *s, int argc, t_atom *argv) {
    int nmax;
    MTXQhull *xo = (MTXQhull *) pd_new (mtx_qhull_class);
    xo->outl = outlet_new (&xo->x_obj, gensym("matrix"));
    xo->zh=0;
    xo->list=0;
    xo->size=0;
    xo->iter=0;
    if (argc>0)
        mTXQhullSetIterations(xo,atom_getfloat(argv));
    return ((void *) xo);
}


static void mTXQhullMatrix(MTXQhull *xo, t_symbol *s,
        int argc, t_atom *argv) {
    int i;
    int rows=atom_getint(argv++);
    int columns=atom_getint(argv++);
    int size=rows*columns;
    int in_size=argc-2;
    float *x;
    float *y;
    float *z;

    /* size check */
    if (!size) 
        post("mtx_qhull: invalid dimensions");
    else if (in_size<size) 
        post("mtx_qhull: sparse matrix not yet supported: use \"mtx_check\"");
    else if ((rows<4)||(columns!=3))
        post("mtx_qhull: requires an L x 3 matrix with at least L>=4");
    else {
        xo->zh = (zhull_t*)malloc(sizeof(zhull_t));
        x=(float*)malloc(rows*sizeof(float));
        y=(float*)malloc(rows*sizeof(float));
        z=(float*)malloc(rows*sizeof(float));
        if ((x!=0)&&(y!=0)&&(z!=0)&&(xo->zh!=0)) {
            for (i=0; i<rows; i++) {
                x[i]=atom_getfloat(argv++);
                y[i]=atom_getfloat(argv++);
                z[i]=atom_getfloat(argv++);
            }
            *(xo->zh)=zhullInitPoints(x,y,z,rows);	
            calculateZHull(xo->zh,xo->iter);
            free(x);
            free(y);
            free(z);
            xo->hull_size=getLength(xo->zh->facets);
            if (xo->list==0) {
                xo->list = (t_atom*)malloc((xo->hull_size*3+2)*sizeof(t_atom));
            } else {
                xo->list = (t_atom*)realloc(xo->list,(xo->hull_size*3+2)*sizeof(t_atom));
            }
            if(xo->list!=0) {
                xo->size=(xo->hull_size*3+2);
                SETFLOAT(xo->list,(float)xo->hull_size);
                SETFLOAT(xo->list+1,(float)3);
                for (i=0; i<xo->hull_size; i++) {
                    SETFLOAT(xo->list+2+3*i, (float)getEntry(
                                getFacetByIndex(xo->zh->facets,i)->corners,0)+1);
                    SETFLOAT(xo->list+3+3*i, (float)getEntry(
                                getFacetByIndex(xo->zh->facets,i)->corners,1)+1);
                    SETFLOAT(xo->list+4+3*i, (float)getEntry(
                                getFacetByIndex(xo->zh->facets,i)->corners,2)+1);
                }
                outlet_anything(xo->outl, gensym("matrix"), 
                        xo->size, xo->list);
                freeZhull(xo->zh);
                free(xo->zh);
                xo->zh=0;
            }
            else {
                post("mtx_qhull: memory problem, no operation!");
                xo->size=0;
                freeZhull(xo->zh);
                free(xo->zh);
                xo->zh=0;
            } 
        }
        else {
            if(x!=0)
                free(x);
            if(y!=0)
                free(y);
            if(z!=0)
                free(z);
            if(xo->zh!=0)
                free(xo->zh);            
            xo->zh=0;
            post("mtx_qhull: memory error, no operation!");
        }

    }
}

void mtx_qhull_setup (void) {
    mtx_qhull_class = class_new (
            gensym("mtx_qhull"),
            (t_newmethod) newMTXQhull,
            (t_method) deleteMTXQhull,
            sizeof(MTXQhull),
            CLASS_DEFAULT, A_GIMME, 0);
    class_addmethod(mtx_qhull_class, (t_method) mTXQhullMatrix, gensym("matrix"), A_GIMME, 0);
    class_addfloat(mtx_qhull_class, (t_method) mTXQhullSetIterations);
}

void iemtx_qhull_setup(void){
    mtx_qhull_setup();
}


