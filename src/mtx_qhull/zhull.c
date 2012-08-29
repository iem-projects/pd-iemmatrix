#include "zhull.h"
/*
 *  zhull
 *
 *  own qhull algorithm implementation
 *
 * Copyright (c) 2012, Franz Zotter
 * with friendly help from 
 * IOhannes zmoelnig
 * IEM, Graz, Austria
 * 
 * own Implementation after the QHULL algorithm
 * that is documented in
 * Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., 
 * "The Quickhull algorithm for convex hulls," ACM Trans.
 * on Mathematical Software, 22(4):469-483, Dec 1996, 
 * http://www.qhull.org
 *
 */





/* facets, memory things */

static void freeFacet(facet_t *facet) {
    /*    printf("deleting facet %li\n",facet);
          printList(facet->corners);
          printList(facet->outsideset);
          printList(facet->insideset);
          printList(facet->neighbors);*/
    freeList(&(facet->corners));
    freeList(&(facet->outsideset));
    freeList(&(facet->insideset));
    freeList(&(facet->neighbors));
}

static list_t appendNewFacets(zhull_t * const zh, const size_t num_facets) {
    facet_t *new_facet;
    index_t i;
    entry_t e0=entry_makeIndex(0);
    list_t new_facets;
    
    new_facets=initConstantList(e0, num_facets);

    for (i=0; i<getLength(new_facets); i++) {
        new_facet = (facet_t*) malloc(sizeof(facet_t));
        if (new_facet==0)
            break;
        new_facet->neighbors=emptyList();
        new_facet->corners=emptyList();
        new_facet->outsideset=emptyList();
        new_facet->insideset=emptyList();
        new_facet->maxdistance=0;
        new_facet->farthest_outside_point=0;
        //        printf("%li, %li\n", new_facet, entry_makePointer(new_facet));
        setEntry(new_facets,i,entry_makePointer(new_facet));
        //        printf("created facet %li\n",new_facet);
    }
    appendListToList(&(zh->facets),new_facets);
    return new_facets;
}

/* facets, interface */
/*entry_t getFacetCornerByIndex(const facet_t *facet,
  index_t index_corner) {
  if (facet!=0)
  return getEntry(facet->corners, 
  index_corner%getLength(facet->corners));
  else
  return 0;
  }
  */

static facet_t *getFacetByIndex(const list_t facets, 
                         const index_t index) {
  entry_t e=getEntry(facets,index);
  return ((facet_t*)entry_getPointer(&e));
}

index_t getTriangleCorner(const zhull_t * const zh,
        const index_t triangle_idx,
        const index_t corner_idx) {
  if (triangle_idx<getLength(zh->facets)) {
    entry_t e=getEntry(getFacetByIndex(zh->facets,triangle_idx)->corners,corner_idx);
    return entry_getIndex(&e);
  } else
    return 0;
}

static facet_t *getFacetByPointer(const entry_t e) {
  return (facet_t*) entry_getPointer(&e);
}

static void getHorizonEdgeByIndex(index_t *corners,
        const list_t horizon_fcts, const list_t horizon_fcts_edges,
        const list_t other_horizon_edges, const index_t index) {
    index_t i=(index+getLength(horizon_fcts_edges))
        %getLength(horizon_fcts_edges);

    entry_t  e = getEntry(horizon_fcts_edges,i);

    facet_t *f = getFacetByIndex(horizon_fcts,i);
    index_t j= entry_getIndex(&e);

    if (f==0) {
      e=getEntry(horizon_fcts_edges,i);
      corners[0]=entry_getIndex(&e);
      e=getEntry(other_horizon_edges,i);
      corners[1]=entry_getIndex(&e);
    } else {
      e=getEntry(f->corners, j);
      corners[0]=entry_getIndex(&e);
      if (getLength(f->corners)>0) {
        e=getEntry(f->corners, (j+1) % getLength(f->corners));
        corners[1]=entry_getIndex(&e);
      } else 
        corners[1]=corners[0];
    }
}


static void removeFacetByPointer(zhull_t * const zh, 
        facet_t * const pointer) {
  removeValueFromList(&(zh->facets), entry_makePointer(pointer));
  removeValueFromList(&(zh->facets_with_outsidepoints), 
                      entry_makePointer(pointer));
  removeValueFromList(&(zh->facets_with_insidepoints), 
                      entry_makePointer(pointer));
  freeFacet(pointer);
} 

static void removeFacetByPointerList(zhull_t * const zh, 
        const list_t pointers) {
    index_t i;
    for (i=0; i<getLength(pointers); i++) {
        removeFacetByPointer(zh,
                getFacetByIndex(pointers,i));
    }
} 

static void removeFacetByIndex(zhull_t * const zh, 
        const index_t index) {
    removeFacetByPointer(zh,
            getFacetByIndex(zh->facets, index));
} 

static void removeFacetByIndexList(zhull_t * const zh, 
        const list_t indices) {
    facet_t *f; 
    index_t i;
    for (i=0; i<getLength(indices); i++) {
      entry_t e = getEntry(indices,i);
      f = getFacetByIndex(zh->facets, entry_getIndex(&e));
      removeFacetByPointer(zh,f);
    }
}

static void freeFacets(zhull_t * const zh) {
    int i;
    facet_t *f;
    if (getLength(zh->facets)>0) {
        for (i=0; i<getLength(zh->facets); i++) {
            f=getFacetByIndex(zh->facets,i);
            //            printf("deleting facet %li\n",i);
            freeFacet(f);
        }
        freeList(&(zh->facets));
    }
}

void freeZhull(zhull_t *zh) {
    if (zh!=0) {
        //        printf("free zhull\n");
        freeFacets(zh);
        freeList(&(zh->facets_with_insidepoints));
        freeList(&(zh->facets_with_outsidepoints));
        freePoints(&(zh->pts));
    }
}


// ***********************************
//
// interface




static line_t getLine(const zhull_t * const zh, 
        const facet_t * const f,
        const index_t corner) {
    vector_t a, b;
    entry_t e = getEntry(f->corners,  (corner)%getLength(f->corners));
    a=getPoint(zh->pts, entry_getIndex(&e));
    e=getEntry(f->corners,(corner+1)%getLength(f->corners));
    b=getPoint(zh->pts, entry_getIndex(&e));
    return lineFromTwoPoints(a,b);
}

zhull_t zhullInitPoints(const float *x, const float *y, 
        const float *z, const size_t num_points) {
    zhull_t zh;
    zh.pts=initPoints(x,y,z,num_points);
    zh.facets=emptyList();
    zh.facets_with_outsidepoints=emptyList();
    zh.facets_with_insidepoints=emptyList();
    return zh;
}



static void dividePointsBetweenNewFacets (
        zhull_t * const zh, const list_t assoc, 
        const list_t new_facets) {
    index_t i,j;
    facet_t *f;
    list_t point_inside_facet_list=emptyList();
    float d;
    entry_t e;
    for (i=0; i<getLength(assoc); i++) {
      e=getEntry(assoc,i);
      int idx = entry_getIndex(&e);
      for (j=0; j<getLength(new_facets); j++) {
        f=getFacetByIndex(new_facets,j);
        d=distancePointPlane(getPoint(zh->pts,idx),
                             f->plane); 
        if (d>=TOL_OUTSIDEPOINT) 
          break;
        else if (d>=TOL_INSIDEPOINT) {
          appendToList(&point_inside_facet_list,entry_makePointer(f));
        }
      }
      if (d>=TOL_OUTSIDEPOINT) {
        appendToList(&(f->outsideset), e);
        if (notInList(entry_makePointer(f),zh->facets_with_outsidepoints)) 
          appendToList(&(zh->facets_with_outsidepoints),entry_makePointer(f));
        if (f->maxdistance<d) {
          f->maxdistance=d;
          f->farthest_outside_point=idx;
        }
      }
      else {
        if (getLength(point_inside_facet_list)>0) {
          for (j=0; j<getLength(point_inside_facet_list); j++) {
            f=getFacetByIndex(point_inside_facet_list,j);
            
            if (notInList(e,f->insideset)) {
              appendToList(&(f->insideset),e);
            }
          }
          appendListToList(&(zh->facets_with_insidepoints),point_inside_facet_list);
          uniquefyListEntries(&(zh->facets_with_insidepoints));
        }
      }
    }
    freeList(&point_inside_facet_list);
}

static void zhullInitialFacets(zhull_t *zh) {
    list_t assoc = emptyList();
    list_t new_facets = emptyList(); 
    index_t i;
    index_t idx[3]={0,1,2};
    list_t list=initListIndex(idx,3);
    if (getNumPoints(zh->pts)>= 3) {
      //      printf("initListFromTo: %d..%d\n", 0,getNumPoints(zh->pts)-1);
      assoc = initListFromTo(0,getNumPoints(zh->pts)-1);
      //     printList(assoc);
        new_facets = appendNewFacets(zh,2);
        if (getLength(new_facets)==2) {
            do {  // circumvent coincident points
              if (distancePointPoint(
                                     getPoint(zh->pts,idx[0]),
                                     getPoint(zh->pts,idx[1])
                                     )>TOL_DEGENERATE)
                    break;
                else
                    idx[1]++;
            } while (idx[1]<getNumPoints(zh->pts));
            if (idx[1]<getNumPoints(zh->pts)) { 
                do { // circumvent degenerate triangles
                    list=initListIndex(idx,3);
                    if (lengthVector(normalOfListedPoints(zh->pts,list))>TOL_DEGENERATE)
                        break;
                    else {
                        idx[2]++; 
                        freeList(&list);
                    }
                } while (idx[2]<getNumPoints(zh->pts));
                if (idx[2]<getNumPoints(zh->pts)) {
                    getFacetByIndex(new_facets,0)->corners = list;
                    appendListToList(&(zh->used_pts),list);
                    list=initListIndex(idx,3);
                    reverseList(&list);
                    getFacetByIndex(new_facets,1)->corners = list;
                    getFacetByIndex(new_facets,0)->neighbors = 
                      initConstantList(entry_makePointer(getFacetByIndex(new_facets,1)),3);
                    getFacetByIndex(new_facets,1)->neighbors = 
                      initConstantList(entry_makePointer(getFacetByIndex(new_facets,0)),3);
                    for (i=0; i<2; i++) {
                        getFacetByIndex(new_facets,i)->plane = 
                            planeFromListedPoints(zh->pts, getFacetByIndex(new_facets,i)->corners);
                        getFacetByIndex(new_facets,i)->outsideset = emptyList();
                        getFacetByIndex(new_facets,i)->insideset = emptyList();
                        getFacetByIndex(new_facets,i)->maxdistance = 0.0f;
                        //printf("removing facests\n");
                        //printList(getFacetByIndex(new_facets,i)->corners);
                        removeValueListFromList(&assoc,getFacetByIndex(new_facets,i)->corners);
                    }
                    //printf("dividePoints...");
                    //printList(assoc);
                    dividePointsBetweenNewFacets(zh, assoc, new_facets);
                }
            }
        }
        freeList(&new_facets);
        freeList(&assoc); 
    }
}

static void printHorizonEdges(list_t *horizon_fcts, 
        list_t *horizon_fcts_edges, 
        list_t *other_horizon_edges) {
    index_t i;
    index_t c1[2];
    list_t list_for_printing=emptyList();

    printf("horizon : ");
    printList(*horizon_fcts);
    printList(*horizon_fcts_edges);
    printList(*other_horizon_edges);
    printf("\n horizon edges: ");
    for (i=0; i<getLength(*horizon_fcts_edges); i++) {
        getHorizonEdgeByIndex(c1, *horizon_fcts, *horizon_fcts_edges,
               *other_horizon_edges, i);
        printf(", %d->%d",c1[0],c1[1]);
    }
    printf("\n");
}


static void sortHorizonEdges(list_t *horizon_fcts, 
        list_t *horizon_fcts_edges, 
        list_t *other_horizon_edges) {
    index_t i,j;
    facet_t *fi;
    entry_t ei;
    index_t c1[2];
    index_t c2[2];

    if (getLength(*horizon_fcts_edges)==0)
        return;
    for (i=0; i<getLength(*horizon_fcts_edges)-1; i++) {
        getHorizonEdgeByIndex(c1, *horizon_fcts, *horizon_fcts_edges,
               *other_horizon_edges, i);
        for (j=i+1; j<getLength(*horizon_fcts_edges); j++) {
            getHorizonEdgeByIndex(c2, *horizon_fcts, *horizon_fcts_edges, 
                    *other_horizon_edges, j);
            if ((c1[1]==c2[0])&&(c1[0]!=c2[1])) { 
                //found edge continuation w/o u-turn: swap positions
                ei=getEntry(*horizon_fcts,j);
                setEntry(*horizon_fcts,j,getEntry(*horizon_fcts,i+1));
                setEntry(*horizon_fcts,i+1,ei);
                ei=getEntry(*horizon_fcts_edges,j);
                setEntry(*horizon_fcts_edges,j,getEntry(*horizon_fcts_edges,i+1));
                setEntry(*horizon_fcts_edges,i+1,ei);
                ei=getEntry(*other_horizon_edges,j);
                setEntry(*other_horizon_edges,j,getEntry(*other_horizon_edges,i+1));
                setEntry(*other_horizon_edges,i+1,ei);
                break;
            }
        }
        if (j==getLength(*horizon_fcts_edges)) {
            // edge continuation requires u-turn 
            for (j=i+1; j<getLength(*horizon_fcts_edges); j++) {
                getHorizonEdgeByIndex(c2, *horizon_fcts, *horizon_fcts_edges, 
                        *other_horizon_edges, j);
                if (c1[1]==c2[0]) { 
                    //found edge continuation w/o return: swap positions
                    ei=getEntry(*horizon_fcts,j);
                    setEntry(*horizon_fcts,j,getEntry(*horizon_fcts,i+1));
                    setEntry(*horizon_fcts,i+1,ei);
                    ei=getEntry(*horizon_fcts_edges,j);
                    setEntry(*horizon_fcts_edges,j,getEntry(*horizon_fcts_edges,i+1));
                    setEntry(*horizon_fcts_edges,i+1,ei);
                    ei=getEntry(*other_horizon_edges,j);
                    setEntry(*other_horizon_edges,j,getEntry(*other_horizon_edges,i+1));
                    setEntry(*other_horizon_edges,i+1,ei);
                    break;
                }
            }
        }
    }
}

static void removeVisibleFacetsGetHorizonAndAvailablePoints(
        zhull_t * const zh, index_t point_index, 
        facet_t *f,
        list_t *horizon_fcts, 
        list_t *horizon_fcts_edges,
        list_t *other_horizon_edges,
        list_t *avail_points) {

    index_t j,k;
    facet_t *n;
    float d;
    list_t visible_fcts = emptyList();
    list_t fcts_to_visit = emptyList();
    list_t fcts_visited = emptyList();
//    list_t list_for_printing = emptyList();

    *avail_points = emptyList();
    *horizon_fcts = emptyList();
    *horizon_fcts_edges = emptyList();
    *other_horizon_edges = emptyList();

    d=distancePointPlane(getPoint(zh->pts,point_index),f->plane);
//    printf("distance %5.2f\n",d);
    appendToList(&fcts_to_visit,entry_makePointer(f));
    if (d>=TOL_OUTSIDEPOINT) { 
        while (getLength(fcts_to_visit)>0) { 
            // visiting only visible facets
            // horizon: edges to invisible or coincident neighbors
            f=getFacetByIndex(fcts_to_visit,0);
            appendToList(&visible_fcts,entry_makePointer(f));
            appendListToList(avail_points, f->outsideset);
            for (j=0; j<getLength(f->neighbors); j++) {
                n=getFacetByIndex(f->neighbors,j);
                d=distancePointPlane(getPoint(zh->pts,point_index),n->plane);
                if (d>=TOL_OUTSIDEPOINT) {  // visit visible neighbors
                    appendToList(&fcts_to_visit,entry_makePointer(n));
                }
                else { // horizon: coincident or invisible neighbors 
                  k=findValueInList(getEntry(f->corners,(j+1)%getLength(f->corners)),
                                    n->corners);
                  appendToList(horizon_fcts,entry_makePointer(n));
                  appendToList(horizon_fcts_edges,entry_makeIndex(k));
                  appendToList(other_horizon_edges,entry_makeIndex(getNumPoints(zh->pts)));
                }
            }
            removeValueFromList(&fcts_to_visit,entry_makePointer(f));
            appendToList(&fcts_visited,entry_makePointer(f));
            removeValueListFromList(&fcts_to_visit,fcts_visited);
        }
//        printf("removing facets\n");
//        list_for_printing=findValueListInList(visible_fcts,zh->facets);
//        printList(list_for_printing);
//        freeList(&list_for_printing);
        removeFacetByPointerList(zh,visible_fcts);
        freeList(&visible_fcts);
        freeList(&fcts_to_visit);
        freeList(&fcts_visited);
        sortHorizonEdges(horizon_fcts, horizon_fcts_edges, other_horizon_edges);
        //printHorizonEdges(horizon_fcts,horizon_fcts_edges,other_horizon_edges);
    }
    else if (d>=TOL_INSIDEPOINT) {
        // all coincident surfaces shall be removed
        // horizon might not be defined by facets
        while (getLength(fcts_to_visit)>0) {
            f=getFacetByIndex(fcts_to_visit,0);
            appendToList(&visible_fcts,entry_makePointer(f));
            appendListToList(avail_points, f->outsideset);
            appendListToList(avail_points, f->insideset);
            for (j=0;j<getLength(f->neighbors);j++) {
                n=getFacetByIndex(f->neighbors,j);
                d=distancePointPlane(getPoint(zh->pts,point_index),n->plane);
                if (d>=TOL_INSIDEPOINT) { // coincident facet 
                    if (notInList(entry_makePointer(n),visible_fcts)) 
                        appendToList(&fcts_to_visit,entry_makePointer(n));
                    if ((innerProduct(f->plane.normal,n->plane.normal)<
                            -1.0f+TOL_DEGENERATE)&&
                            (distancePointLineOnPlane(getPoint(zh->pts,point_index),
                                                      getLine(zh,f,j),
                                                      f->plane)<TOL_INSIDEPOINT)) {
                        // coincident facets with opposite surface orientation
                        // yield an edge to keep despite all facets will be removed 
                        // as soon as edge is invisible to point
                      appendToList(horizon_fcts,entry_makePointer(0));
                      appendToList(other_horizon_edges,
                                   getEntry(f->corners,j));
                      appendToList(horizon_fcts_edges,
                                   getEntry(f->corners,(j+1)%getLength(f->corners)));                
                    }
                }
                else { // invisible facet forms horizon that persists
                    k=findValueInList(getEntry(f->corners,(j+1)%getLength(f->corners)),
                                      n->corners);

                    appendToList(horizon_fcts,entry_makePointer(n));
                    appendToList(horizon_fcts_edges,entry_makeIndex(k));
                    appendToList(other_horizon_edges,entry_makeIndex(getNumPoints(zh->pts)));
                }
            }
            removeValueFromList(&fcts_to_visit,entry_makePointer(f));
            appendToList(&fcts_visited,entry_makePointer(f));
            removeValueListFromList(&fcts_to_visit,fcts_visited);
        }
//        printf("removing facets\n");
//        list_for_printing=findValueListInList(visible_fcts,zh->facets);
//        printList(list_for_printing);
//        freeList(&list_for_printing);
        removeFacetByPointerList(zh,visible_fcts);
        sortHorizonEdges(horizon_fcts, horizon_fcts_edges,other_horizon_edges);
        //printHorizonEdges(horizon_fcts,horizon_fcts_edges,other_horizon_edges);
        freeList(&visible_fcts);
        freeList(&fcts_to_visit);
        freeList(&fcts_visited);
    }
}

static void initNewFacets(zhull_t *zh, index_t point_index,
        list_t new_facets, list_t horizon_fcts,
        list_t horizon_fcts_edges, list_t other_horizon_edges) {
    index_t i,j;
    entry_t entry_array[3];
    index_t array[3];
    index_t temp;
    entry_t e;
    facet_t*f;

    //array[0]=entry_makeIndex(point_index);
    for (i=0; i<getLength(new_facets); i++) {
        // corners
        getHorizonEdgeByIndex(array,
                              horizon_fcts, 
                              horizon_fcts_edges,
                              other_horizon_edges,
                              i);
        temp=array[1];
        array[1]=array[0];
        array[0]=temp;

        array[2]=point_index;
        f=getFacetByIndex(new_facets,i);
        //        printf("facets=%p\n", f);
        f->corners=initListIndex(array,3);

        // neighbors
        // previous new neighbor
        j=(getLength(horizon_fcts)+i-1) % getLength(horizon_fcts);
        entry_array[1]=getEntry(new_facets,j);
        // next new neighbor
        j=(i+1) % getLength(horizon_fcts);
        entry_array[2]=getEntry(new_facets,j);
        // old neighbor
        e=getEntry(horizon_fcts,i);
        if (entry_getPointer(&e)!=0) {
          e=getEntry(horizon_fcts_edges,i);
          entry_array[0]=getEntry(horizon_fcts,i);
          setEntry(getFacetByIndex(horizon_fcts,i)->neighbors,
                   entry_getIndex(&e),
                   getEntry(new_facets,i));
        } else {
            // registring at new neighbor where there
            // is no old one: degenerate 2D case
           for (j=0;j<getLength(horizon_fcts_edges);j++) {
             entry_t e1=getEntry(horizon_fcts_edges,i);
             entry_t e2=getEntry(horizon_fcts_edges,j);
             entry_t e3=getEntry(other_horizon_edges,i);
             entry_t e4=getEntry(other_horizon_edges,j);



             if (entry_equals(&e2, &e3)&&entry_equals(&e1,&e4))
               break;
           }
           entry_array[0]=getEntry(new_facets,j);
           facet_t*fp= getFacetByIndex(new_facets, j);
           list_t neighbors =getFacetByIndex(new_facets,j)->neighbors;
           setEntry(neighbors,
                    0, getEntry(new_facets,i));
        }

        getFacetByIndex(new_facets,i)->neighbors=
            initList(entry_array,3);

        // removing inside points at (potential) horizon facet if point index is one
        if (getFacetByIndex(horizon_fcts,i)!=0) {
          removeValueFromList(&(getFacetByIndex(horizon_fcts,i)->insideset),entry_makeIndex(point_index));
            if (getLength(getFacetByIndex(horizon_fcts,i)->insideset)==0) {
                removeValueFromList(&(zh->facets_with_insidepoints),
                                    entry_makePointer(getFacetByIndex(horizon_fcts,i)));
            }
        }

        // initializing normal vectors and lists
        getFacetByIndex(new_facets,i)->plane = 
            planeFromListedPoints(zh->pts, 
                    getFacetByIndex(new_facets,i)->corners);

//        printf("new facet %d ",findValueInList(getEntry(new_facets,i),zh->facets));
//        printFacet(zh, getFacetByIndex(new_facets,i));
    }
}

static void makePyramidFacetsToHorizon(zhull_t *zh, index_t point_index, 
        list_t horizon_fcts, list_t horizon_fcts_edges,
        list_t other_horizon_edges, list_t avail_points) {
    list_t new_facets = appendNewFacets(zh, getLength(horizon_fcts_edges));
    //    printf("making new pyramid of %d facets\n",getLength(horizon_fcts_edges));
    initNewFacets(zh,point_index,new_facets,horizon_fcts,horizon_fcts_edges,other_horizon_edges);
    appendToList(&(zh->used_pts),entry_makeIndex(point_index));
/*        printf("available points: ");
        printList(avail_points);
        printf("new facets : ");
        printList(new_facets);*/
    dividePointsBetweenNewFacets(zh, avail_points, new_facets);
    freeList(&new_facets);
}

static appendExteriorPoints(zhull_t *zh) {
    index_t i;
    vector_t center = initVector(0.0f,0.0f,0.0f);
    list_t facet_delete_list=emptyList();
    facet_t *f;
    center=averageListedPoints(zh->pts,zh->used_pts);
    printf("central point\n");
    printVector(center);
    printf("\n");
    for (i=0; i<getLength(zh->facets); i++) {
        f=getFacetByIndex(zh->facets,i);
        printf("distance of plane %d, d=%5.2f\n",i,
                distancePointPlane(center,f->plane));
        if (distancePointPlane(center,f->plane)>-0.5f) {
            appendToList(&facet_delete_list,entry_makePointer(f));
        }
    }
    printList(facet_delete_list);
    removeFacetByPointerList(zh,facet_delete_list);
    freeList(&facet_delete_list);
}

int calculateZHull(zhull_t *zh) {
    index_t fli=0;
    index_t pi;
    facet_t *f;
    list_t outsideset;
    int cnt=0;
    int maxit=getNumPoints(zh->pts);
    list_t horizon_fcts=emptyList();
    list_t horizon_fcts_edges=emptyList();
    list_t other_horizon_edges=emptyList();
    list_t available_points=emptyList();
    entry_t e;


//    if (maxit>MAXIT)
//        maxit=MAXIT;
    if (getNumPoints(zh->pts)!=0){
        zhullInitialFacets(zh);
        //printZhull(zh);
        while(((getLength(zh->facets_with_insidepoints)>0)
                    ||(getLength(zh->facets_with_outsidepoints)>0))
                &&(cnt++<maxit))  {
          //            printf("//////////////// ITERATION %d ///////\n",cnt);
            if (getLength(zh->facets_with_insidepoints)>0) {
                fli%=getLength(zh->facets_with_insidepoints);
                f=getFacetByIndex(zh->facets_with_insidepoints,fli);
                e=getEntry(f->insideset, 0);
                pi=entry_getIndex(&e);
//                printf("insidepoint\n");
//                printList(zh->facets_with_insidepoints);
            }
            else {
                fli%=getLength(zh->facets_with_outsidepoints);
                f=getFacetByIndex(zh->facets_with_outsidepoints,fli);
                pi=f->farthest_outside_point;
            }
//            printf("point %d\n",pi);
            removeVisibleFacetsGetHorizonAndAvailablePoints(zh,pi,f,
                    &horizon_fcts, &horizon_fcts_edges,&other_horizon_edges,
                    &available_points);
            removeValueFromList(&available_points, entry_makeIndex(pi));
            makePyramidFacetsToHorizon(zh,pi,horizon_fcts,horizon_fcts_edges,
                    other_horizon_edges,available_points);
            //            printZhull(zh);
            freeList(&horizon_fcts);
            freeList(&horizon_fcts_edges);
            freeList(&other_horizon_edges);
            freeList(&available_points);
            fli++;
        }
//        appendExteriorPoints(zh);
    }
    return cnt;
}

void printZhull(const zhull_t * const zh) {
    index_t fi;
    list_t indices = emptyList();
/*     printf("zhull from %d points\n", getNumPoints(zh->pts));
    printf("facets with outsidepoints: ");
    indices=findValueListInList(zh->facets_with_outsidepoints,zh->facets);
    printList(indices);
    freeList(&indices);
    printf("facets with insidepoints: ");
    indices=findValueListInList(zh->facets_with_insidepoints,zh->facets);
    printList(indices);
    freeList(&indices);
    */
    printf("zhull has %d facets\n", getLength(zh->facets));
    for (fi=0; fi<getLength(zh->facets); fi++) {
        printf("facet %d<%d>: ",fi,getFacetByIndex(zh->facets,fi));
        printFacet(zh,getFacetByIndex(zh->facets,fi));
    }
}

void printFacet(const zhull_t * const zh, 
        const facet_t * const f) {
    list_t indices=emptyList();
    indices=findValueListInList(f->neighbors,zh->facets);
    printf("plane: ");
    printPlane(f->plane);
    printf("\n");
    printf("corners: ");
    printList(f->corners);
    printf("outsideset: ");
    printList(f->outsideset);
    printf("insideset: ");
    printList(f->insideset);
    printf("neighbors: ");
    printList(indices);
    freeList(&indices);
    printf("pt %d with maxdist %5.2f\n",f->farthest_outside_point, f->maxdistance);
}

