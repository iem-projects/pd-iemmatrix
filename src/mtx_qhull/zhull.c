#include "zhull.h"
/* facets, memory things */

void freeFacet(facet_t *facet) {
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

list_t appendNewFacets(zhull_t * const zh, const size_t num_facets) {
    facet_t *new_facet;
    list_t new_facets=initConstantList(0,num_facets);
    index_t i;
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
        //        printf("%li, %li\n", new_facet, (entry_t)new_facet);
        setEntry(new_facets,i,(entry_t)new_facet);
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

facet_t *getFacetByIndex(const list_t facets, 
        const index_t index) {
    return (facet_t*) getEntry(facets,index);
}

index_t getTriangleCorner(const zhull_t * const zh,
        const index_t triangle_idx,
        const index_t corner_idx) {
    if (triangle_idx<getLength(zh->facets)) 
        return getEntry(getFacetByIndex(zh->facets,triangle_idx)->corners,corner_idx);
    else
        return 0;
}

facet_t *getFacetByPointer(const entry_t ptr) {
    return (facet_t*) ptr;
}

void getHorizonEdgeByIndex(index_t *corners,
        const list_t horizon_fcts, const list_t horizon_fcts_edges,
        const list_t other_horizon_edges, const index_t index) {
    index_t i=(index+getLength(horizon_fcts_edges))
        %getLength(horizon_fcts_edges);
    facet_t *f = getFacetByIndex(horizon_fcts,i);
    index_t j=getEntry(horizon_fcts_edges,i);

    if (f==0) {
        corners[0]=getEntry(horizon_fcts_edges,i);
        corners[1]=getEntry(other_horizon_edges,i);
    } else {
        corners[0]=getEntry(f->corners, j);
        if (getLength(f->corners)>0) 
            corners[1]=getEntry(f->corners, (j+1) % getLength(f->corners));
        else 
            corners[1]=corners[0];
    }
}


void removeFacetByPointer(zhull_t * const zh, 
        facet_t * const pointer) {
    removeValueFromList(&(zh->facets), (entry_t)pointer);
    removeValueFromList(&(zh->facets_with_outsidepoints), 
            (entry_t)pointer);
    removeValueFromList(&(zh->facets_with_insidepoints), 
            (entry_t)pointer);
    freeFacet((facet_t*)pointer);
} 

void removeFacetByPointerList(zhull_t * const zh, 
        const list_t pointers) {
    index_t i;
    for (i=0; i<getLength(pointers); i++) {
        removeFacetByPointer(zh,
                getFacetByIndex(pointers,i));
    }
} 

void removeFacetByIndex(zhull_t * const zh, 
        const index_t index) {
    removeFacetByPointer(zh,
            getFacetByIndex(zh->facets, index));
} 

void removeFacetByIndexList(zhull_t * const zh, 
        const list_t indices) {
    facet_t *f; 
    index_t i;
    for (i=0; i<getLength(indices); i++) {
        f = getFacetByIndex(zh->facets, getEntry(indices,i));
        removeFacetByPointer(zh,f);
    }
} 


void freeFacets(zhull_t * const zh) {
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




line_t getLine(const zhull_t * const zh, 
        const facet_t * const f,
        const index_t corner) {
    vector_t a, b;
    a=getPoint(zh->pts,getEntry(f->corners,  (corner)%getLength(f->corners)));
    b=getPoint(zh->pts,getEntry(f->corners,(corner+1)%getLength(f->corners)));
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



void dividePointsBetweenNewFacets (
        zhull_t * const zh, const list_t assoc, 
        const list_t new_facets) {
    index_t i,j;
    facet_t *f;
    list_t point_inside_facet_list=emptyList();
    float d;
    for (i=0; i<getLength(assoc); i++) {
        for (j=0; j<getLength(new_facets); j++) {
            f=getFacetByIndex(new_facets,j);
            d=distancePointPlane(getPoint(zh->pts,getEntry(assoc,i)),
                    f->plane); 
            if (d>=TOL_OUTSIDEPOINT) 
                break;
            else if (d>=TOL_INSIDEPOINT) {
                appendToList(&point_inside_facet_list,(entry_t)f);
            }
        }
        if (d>=TOL_OUTSIDEPOINT) {
            appendToList(&(f->outsideset), getEntry(assoc,i));
            if (notInList((entry_t)f,zh->facets_with_outsidepoints)) 
                appendToList(&(zh->facets_with_outsidepoints),(entry_t)f);
            if (f->maxdistance<d) {
                f->maxdistance=d;
                f->farthest_outside_point=getEntry(assoc,i);
            }
        }
        else {
            if (getLength(point_inside_facet_list)>0) {
                for (j=0; j<getLength(point_inside_facet_list); j++) {
                    f=getFacetByIndex(point_inside_facet_list,j);
                    if (notInList(getEntry(assoc,i),f->insideset))
                        appendToList(&(f->insideset),getEntry(assoc,i));
                }
                appendListToList(&(zh->facets_with_insidepoints),point_inside_facet_list);
                uniquefyListEntries(&(zh->facets_with_insidepoints));
            }
        }
    }
    freeList(&point_inside_facet_list);
}

void zhullInitialFacets(zhull_t *zh) {
    list_t assoc = emptyList();
    list_t new_facets = emptyList(); 
    index_t i;
    entry_t idx[3]={0,1,2};
    list_t list=initList(idx,3);
    if (getNumPoints(zh->pts)>= 3) {
        assoc = initListFromTo(0,getNumPoints(zh->pts)-1);
        new_facets = appendNewFacets(zh,2);
        if (getLength(new_facets)==2) {
            do {  // circumvent coincident points
                if (distancePointPoint(getPoint(zh->pts,idx[0]),getPoint(zh->pts,idx[1]))>TOL_DEGENERATE)
                    break;
                else
                    idx[1]++;
            } while (idx[1]<getNumPoints(zh->pts));
            if (idx[1]<getNumPoints(zh->pts)) { 
                do { // circumvent degenerate triangles
                    list=initList(idx,3);
                    if (lengthVector(normalOfListedPoints(zh->pts,list))>TOL_DEGENERATE)
                        break;
                    else {
                        idx[2]++; 
                        freeList(&list);
                    }
                } while (idx[2]<getNumPoints(zh->pts));
                if (idx[2]<getNumPoints(zh->pts)) {
                    getFacetByIndex(new_facets,0)->corners = list;
                    list=initList(idx,3);
                    reverseList(&list);
                    getFacetByIndex(new_facets,1)->corners = list;
                    getFacetByIndex(new_facets,0)->neighbors = 
                        initConstantList((entry_t)getFacetByIndex(new_facets,1),3);
                    getFacetByIndex(new_facets,1)->neighbors = 
                        initConstantList((entry_t)getFacetByIndex(new_facets,0),3);
                    for (i=0; i<2; i++) {
                        getFacetByIndex(new_facets,i)->plane = 
                            planeFromListedPoints(zh->pts, getFacetByIndex(new_facets,i)->corners);
                        getFacetByIndex(new_facets,i)->outsideset = emptyList();
                        getFacetByIndex(new_facets,i)->insideset = emptyList();
                        getFacetByIndex(new_facets,i)->maxdistance = 0.0f;
                        removeValueListFromList(&assoc,getFacetByIndex(new_facets,i)->corners);
                    }
                    dividePointsBetweenNewFacets(zh, assoc, new_facets);
                }
            }
        }
        freeList(&new_facets);
        freeList(&assoc); 
    }
}

void printHorizonEdges(list_t *horizon_fcts, 
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


void sortHorizonEdges(list_t *horizon_fcts, 
        list_t *horizon_fcts_edges, 
        list_t *other_horizon_edges) {
    index_t i,j;
    facet_t *fi;
    index_t ei;
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

void removeVisibleFacetsGetHorizonAndAvailablePoints(
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
    appendToList(&fcts_to_visit,(entry_t)f);
    if (d>=TOL_OUTSIDEPOINT) { 
        while (getLength(fcts_to_visit)>0) { 
            // visiting only visible facets
            // horizon: edges to invisible or coincident neighbors
            f=getFacetByIndex(fcts_to_visit,0);
            appendToList(&visible_fcts,(entry_t)f);
            appendListToList(avail_points, f->outsideset);
            for (j=0; j<getLength(f->neighbors); j++) {
                n=getFacetByIndex(f->neighbors,j);
                d=distancePointPlane(getPoint(zh->pts,point_index),n->plane);
                if (d>=TOL_OUTSIDEPOINT) {  // visit visible neighbors
                    appendToList(&fcts_to_visit,(entry_t)n);
                }
                else { // horizon: coincident or invisible neighbors 
                    k=getEntry(f->corners,(j+1)%getLength(f->corners));
                    k=findValueInList(k,n->corners);
                    appendToList(horizon_fcts,(entry_t)n);
                    appendToList(horizon_fcts_edges,k);
                    appendToList(other_horizon_edges,getNumPoints(zh->pts));
                }
            }
            removeValueFromList(&fcts_to_visit,(entry_t)f);
            appendToList(&fcts_visited,(entry_t)f);
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
//        printHorizonEdges(horizon_fcts,horizon_fcts_edges,other_horizon_edges);
    }
    else if (d>=TOL_INSIDEPOINT) {
        // all coincident surfaces shall be removed
        // horizon might not be defined by facets
        while (getLength(fcts_to_visit)>0) {
            f=getFacetByIndex(fcts_to_visit,0);
            appendToList(&visible_fcts,(entry_t)f);
            appendListToList(avail_points, f->outsideset);
            appendListToList(avail_points, f->insideset);
            for (j=0;j<getLength(f->neighbors);j++) {
                n=getFacetByIndex(f->neighbors,j);
                d=distancePointPlane(getPoint(zh->pts,point_index),n->plane);
                if (d>=TOL_INSIDEPOINT) { // coincident facet 
                    if (notInList((entry_t)n,visible_fcts)) 
                        appendToList(&fcts_to_visit,(entry_t)n);
                    if ((innerProduct(f->plane.normal,n->plane.normal)<
                            -1.0f+TOL_DEGENERATE)&&
                            (distancePointLineOnPlane(getPoint(zh->pts,point_index),
                                                      getLine(zh,f,j),
                                                      f->plane)<TOL_INSIDEPOINT)) {
                        // coincident facets with opposite surface orientation
                        // yield an edge to keep despite all facets will be removed 
                        // as soon as edge is invisible to point
                        appendToList(horizon_fcts,0);
                        appendToList(other_horizon_edges,
                                getEntry(f->corners,j));
                        appendToList(horizon_fcts_edges,
                                getEntry(f->corners,(j+1)%getLength(f->corners)));                
                    }
                }
                else { // invisible facet forms horizon that persists
                    k=getEntry(f->corners,(j+1)%getLength(f->corners));
                    k=findValueInList(k,n->corners);
                    appendToList(horizon_fcts,(entry_t)n);
                    appendToList(horizon_fcts_edges,k);
                    appendToList(other_horizon_edges,getNumPoints(zh->pts));
                }
            }
            removeValueFromList(&fcts_to_visit,(entry_t)f);
            appendToList(&fcts_visited,(entry_t)f);
            removeValueListFromList(&fcts_to_visit,fcts_visited);
        }
//        printf("removing facets\n");
//        list_for_printing=findValueListInList(visible_fcts,zh->facets);
//        printList(list_for_printing);
//        freeList(&list_for_printing);
        removeFacetByPointerList(zh,visible_fcts);
        sortHorizonEdges(horizon_fcts, horizon_fcts_edges,other_horizon_edges);
//        printHorizonEdges(horizon_fcts,horizon_fcts_edges,other_horizon_edges);
        freeList(&visible_fcts);
        freeList(&fcts_to_visit);
        freeList(&fcts_visited);
    }
}

void initNewFacets(zhull_t *zh, index_t point_index,
        list_t new_facets, list_t horizon_fcts,
        list_t horizon_fcts_edges, list_t other_horizon_edges) {
    index_t i,j;
    entry_t array[3];

    array[0]=point_index;
    for (i=0; i<getLength(new_facets); i++) {
        // corners
        getHorizonEdgeByIndex(array,horizon_fcts, horizon_fcts_edges,
                other_horizon_edges,i);
        array[2]=array[1];
        array[1]=array[0];
        array[0]=array[2];
        array[2]=point_index;
        getFacetByIndex(new_facets,i)->corners=initList(array,3);

        // neighbors
        // previous new neighbor
        j=(getLength(horizon_fcts)+i-1) % getLength(horizon_fcts);
        array[1]=getEntry(new_facets,j);
        // next new neighbor
        j=(i+1) % getLength(horizon_fcts);
        array[2]=getEntry(new_facets,j);
        // old neighbor
        if (getEntry(horizon_fcts,i)!=0) {
            array[0]=getEntry(horizon_fcts,i);
            setEntry(getFacetByIndex(horizon_fcts,i)->neighbors,
                    getEntry(horizon_fcts_edges,i), getEntry(new_facets,i));
        }
        else {
            // registring at new neighbor where there
            // is no old one: degenerate 2D case
           for (j=0;j<getLength(horizon_fcts_edges);j++) {
               if ((getEntry(horizon_fcts_edges,i)==getEntry(other_horizon_edges,j))&&
                       (getEntry(horizon_fcts_edges,j)==getEntry(other_horizon_edges,i)))
                   break;
           }
           array[0]=getEntry(new_facets,j);
           setEntry(getFacetByIndex(new_facets,j)->neighbors,
                    0, getEntry(new_facets,i));
        }

        getFacetByIndex(new_facets,i)->neighbors=
            initList(array,3);

        // removing inside points at (potential) horizon facet if point index is one
        if (getFacetByIndex(horizon_fcts,i)!=0) {
            removeValueFromList(&(getFacetByIndex(horizon_fcts,i)->insideset),point_index);
            if (getLength(getFacetByIndex(horizon_fcts,i)->insideset)==0) {
                removeValueFromList(&(zh->facets_with_insidepoints),(entry_t)getFacetByIndex(horizon_fcts,i));
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

void makePyramidFacetsToHorizon(zhull_t *zh, index_t point_index, 
        list_t horizon_fcts, list_t horizon_fcts_edges,
        list_t other_horizon_edges, list_t avail_points) {
    list_t new_facets = appendNewFacets(zh, getLength(horizon_fcts_edges));
//    printf("making new pyramid of %d facets\n",getLength(horizon_fcts_edges));
    initNewFacets(zh,point_index,new_facets,horizon_fcts,horizon_fcts_edges,other_horizon_edges);
/*        printf("available points: ");
        printList(avail_points);
        printf("new facets : ");
        printList(new_facets);*/
    dividePointsBetweenNewFacets(zh, avail_points, new_facets);
    freeList(&new_facets);
}

void calculateZHull(zhull_t *zh,int maxit) {
    index_t fli=0;
    index_t pi;
    facet_t *f;
    list_t outsideset;
    int cnt=0;
    list_t horizon_fcts=emptyList();
    list_t horizon_fcts_edges=emptyList();
    list_t other_horizon_edges=emptyList();
    list_t available_points=emptyList();
    if (maxit>MAXIT)
        maxit=MAXIT;
    if (getNumPoints(zh->pts)!=0){
        zhullInitialFacets(zh);
//        printZhull(zh);
        while(((getLength(zh->facets_with_insidepoints)>0)
                    ||(getLength(zh->facets_with_outsidepoints)>0))
                &&(cnt++<maxit))  {
//            printf("//////////////// ITERATION %d ///////\n",cnt);
            if (getLength(zh->facets_with_insidepoints)>0) {
                fli%=getLength(zh->facets_with_insidepoints);
                f=getFacetByIndex(zh->facets_with_insidepoints,fli);
                pi=getEntry(f->insideset,0);
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
            removeValueFromList(&available_points, pi);
            makePyramidFacetsToHorizon(zh,pi,horizon_fcts,horizon_fcts_edges,
                    other_horizon_edges,available_points);
//            printZhull(zh);
            freeList(&horizon_fcts);
            freeList(&horizon_fcts_edges);
            freeList(&other_horizon_edges);
            freeList(&available_points);
            fli++;
        }
    }
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

