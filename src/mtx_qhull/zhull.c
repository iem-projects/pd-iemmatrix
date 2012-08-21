#define TOL 1e-7
#define MAXIT 1000000


vector_t averageListedPoints(const points_t points, const list_t list) {
    index_t i;
    vector_t m = initVector(0.0f, 0.0f, 0.0f);
    for (i=0; i<getLength(list); i++) {
        m=addVectors(points.v[getEntry(list,i)],m);
    }
    m=scaleVector(m,1.0f/((float)getLength(list)));
    return m;
}

vector_t normalOfListedPoints(const points_t points, const list_t list) {
    vector_t m = averageListedPoints(points,list);
    vector_t n = initVector(0.0f, 0.0f, 0.0f);
    vector_t d1,d2,c;
    index_t i;
    for (i=1; i<getLength(list); i++) {
        d1=subtractVectors(points.v[getEntry(list,i-1)],m);
        d2=subtractVectors(points.v[getEntry(list,i)],m);
        c=crossProduct(d1,d2);
        n=addVectors(c,n);
    }
    return n;
}

plane_t planeFromListedPoints(const points_t points, const list_t list) {
    vector_t p=averageListedPoints(points,list);
    vector_t n=normalOfListedPoints(points,list);
    plane_t pl=initPlane(n,p);
    return pl;
}

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

facet_t *getFacetByPointer(const entry_t ptr) {
    return (facet_t*) ptr;
}

void getHorizonEdgeByIndex(index_t *corners,
        const list_t horizon_fcts, const list_t horizon_fcts_edges, 
        const index_t index) {
    index_t i=(index+getLength(horizon_fcts_edges))
        %getLength(horizon_fcts_edges);
    facet_t *f = getFacetByIndex(horizon_fcts,i);
    index_t j=getEntry(horizon_fcts_edges,i);

    corners[0]=getEntry(f->corners, j);
    if (getLength(f->corners)>0) 
        corners[1]=getEntry(f->corners, (j+1) % getLength(f->corners));
    else 
        corners[1]=corners[0];
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


vector_t getPoint(const zhull_t * const zh, 
        const index_t index) {
    return (vector_t) zh->pts.v[index];
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
        const list_t new_facets, const list_t horizon_fcts) {
    index_t i,j;
    list_t new_and_horizon_fcts=mergeLists(new_facets,horizon_fcts);
    printList(new_and_horizon_fcts);
    uniquefyListEntries(&new_and_horizon_fcts);
    printList(new_and_horizon_fcts);
    float d;
    for (i=0; i<getLength(assoc); i++) {
        printf("p %d\n", getEntry(assoc,i));
        for (j=0; j<getLength(new_and_horizon_fcts); j++) {
            d=distancePointPlane(getPoint(zh,getEntry(assoc,i)),
                    getFacetByIndex(new_and_horizon_fcts,j)->plane); 
 //           printf("distance %5.2f to %d\n",d,findValueInList((entry_t)getFacetByIndex(new_and_horizon_fcts,j),zh->facets));
            /*
               if ((d>=-TOL) && (d<=TOL)) {
               appendToList(&(getFacetByIndex(new_facets,j)->insideset),
               getEntry(assoc,i)); 
               if (notInList(getEntry(new_facets,j),zh->facets_with_insidepoints)) 
               appendToList(&(zh->facets_with_insidepoints),
               getEntry(new_facets,j));
               break;
               }
               else */
            if (d>=-TOL) {
                appendToList(&(getFacetByIndex(new_and_horizon_fcts,j)->outsideset),
                        getEntry(assoc,i));
   //             printf("appended ");
                if (notInList(getEntry(new_and_horizon_fcts,j),zh->facets_with_outsidepoints)) 
                    appendToList(&(zh->facets_with_outsidepoints),getEntry(new_and_horizon_fcts,j));
                if (getFacetByIndex(new_and_horizon_fcts,j)->maxdistance<d) {
                    getFacetByIndex(new_and_horizon_fcts,j)->maxdistance=d;
                    getFacetByIndex(new_and_horizon_fcts,j)->farthest_outside_point=getEntry(assoc,i);
                }
                break;
            }
        }
    }
}

void zhullInitialFacets(zhull_t *zh) {
    list_t assoc = emptyList();
    list_t new_facets = emptyList(); 
    index_t i;
    if (zh->pts.num_points >= 3) {
        assoc = initListFromTo(0,zh->pts.num_points-1);
        new_facets = appendNewFacets(zh,2);
        if (getLength(new_facets)==2) {
            getFacetByIndex(new_facets,0)->corners = initListFromTo(0,2);
            getFacetByIndex(new_facets,1)->corners = initListFromTo(2,0);
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
            dividePointsBetweenNewFacets(zh, assoc, new_facets, emptyList());
        }
        freeList(&new_facets);
        freeList(&assoc); 
    }
}

void printHorizonEdges(list_t *horizon_fcts, 
        list_t *horizon_fcts_edges) {
    index_t i;
    index_t c1[2];

    printf("horizon edges: ");
    for (i=0; i<getLength(*horizon_fcts_edges); i++) {
        getHorizonEdgeByIndex(c1, *horizon_fcts, *horizon_fcts_edges, i);
        printf(", %d->%d",c1[0],c1[1]);
    }
    printf("\n");
}


void sortHorizonEdges(list_t *horizon_fcts, 
        list_t *horizon_fcts_edges) {
    index_t i,j;
    facet_t *fi;
    index_t ei;
    index_t c1[2];
    index_t c2[2];

    if (getLength(*horizon_fcts_edges)==0)
        return;
    for (i=0; i<getLength(*horizon_fcts_edges)-1; i++) {
        getHorizonEdgeByIndex(c1, *horizon_fcts, *horizon_fcts_edges, i);
        for (j=i+1; j<getLength(*horizon_fcts_edges); j++) {
            getHorizonEdgeByIndex(c2, *horizon_fcts, *horizon_fcts_edges, j);
            if (c1[1]==c2[0]) { //found edge continuation: swap positions
                ei=getEntry(*horizon_fcts,j);
                setEntry(*horizon_fcts,j,getEntry(*horizon_fcts,i+1));
                setEntry(*horizon_fcts,i+1,ei);
                ei=getEntry(*horizon_fcts_edges,j);
                setEntry(*horizon_fcts_edges,j,getEntry(*horizon_fcts_edges,i+1));
                setEntry(*horizon_fcts_edges,i+1,ei);
                break;
            }
            // HIER MUSS MEHR LOGIK HINEIN: HORIZONT KANN
            // FEHLER BEINHALTEN
        }
    }
}

void removeVisibleFacetsGetHorizonAndAvailablePoints(
        zhull_t * const zh, index_t point_index, 
        facet_t *facet,
        list_t *horizon_fcts, 
        list_t *horizon_fcts_edges,
        list_t *avail_points) {

    index_t i,j,k;
    facet_t *f, *n;
    float d;
    list_t visible_fcts = emptyList();
    list_t indices_for_printing = emptyList();

    *avail_points = emptyList();
    *horizon_fcts = emptyList();
    *horizon_fcts_edges = emptyList();

    appendToList(&visible_fcts, (entry_t)facet);
    for(i=0;i<getLength(visible_fcts);i++) {
        f=getFacetByIndex(visible_fcts,i);
        appendListToList(avail_points, f->outsideset);
        appendListToList(avail_points, f->insideset);
        for (j=0; j<getLength(f->neighbors); j++) {
            n=getFacetByIndex(f->neighbors,j);
            if (notInList((entry_t) n, visible_fcts)) {
                d=distancePointPlane(getPoint(zh,point_index), 
                        n->plane);
                if (d>-TOL) {
                    if (notInList((entry_t)n,visible_fcts)) 
                        appendToList(&visible_fcts,(entry_t)n);
                }
                else {
                    appendToList(horizon_fcts,(entry_t)n);
                    k=getEntry(f->corners,(j+1)%getLength(f->corners));
                    //		printf("searching for corner with %d in\n",k);
                    //		printList(n->corners);
                    k=findValueInList(k,n->corners);
                    //		printf("found %d\n",k);
                    appendToList(horizon_fcts_edges,k);
                }
            }
        }
    }
    printf("removing facets");
    indices_for_printing=findValueListInList(visible_fcts, zh->facets);
    printList(indices_for_printing);
    removeFacetByPointerList(zh,visible_fcts);
    freeList(&visible_fcts);
    //    printHorizonEdges(horizon_fcts, horizon_fcts_edges);
    sortHorizonEdges(horizon_fcts, horizon_fcts_edges);
    printHorizonEdges(horizon_fcts, horizon_fcts_edges);
}


void initNewFacets(zhull_t *zh, int point_index,
        list_t new_facets, list_t horizon_fcts,
        list_t horizon_fcts_edges) {
    index_t i,j;
    entry_t array[3];

    array[0]=point_index;
    for (i=0; i<getLength(new_facets); i++) {
        //        printf("adding facet %d\n",getEntry(new_facets,i));

        // corners
        getHorizonEdgeByIndex(array,horizon_fcts, horizon_fcts_edges,i);
        array[2]=array[1];
        array[1]=array[0];
        array[0]=array[2];
        array[2]=point_index;
        getFacetByIndex(new_facets,i)->corners=initList(array,3);
        //        printList(getFacetByIndex(new_facets,i)->corners);

        // neighbors
        // previous new neighbor
        j=(getLength(horizon_fcts)+i-1) % getLength(horizon_fcts);
        array[1]=getEntry(new_facets,j);
        // next new neighbor
        j=(i+1) % getLength(horizon_fcts);
        array[2]=getEntry(new_facets,j);
        // old neighbor
        array[0]=getEntry(horizon_fcts,i);
        getFacetByIndex(new_facets,i)->neighbors=
            initList(array,3);

        // registering at old neighbor:
        setEntry(getFacetByIndex(horizon_fcts,i)->neighbors,
                getEntry(horizon_fcts_edges,i), getEntry(new_facets,i));

        // initializing normal vectors and lists
        getFacetByIndex(new_facets,i)->plane = 
            planeFromListedPoints(zh->pts, 
                    getFacetByIndex(new_facets,i)->corners);
    }
}

void makePyramidFacetsToHorizon(zhull_t *zh, index_t point_index, 
        list_t horizon_fcts, list_t horizon_fcts_edges,
        list_t avail_points) {
    list_t new_facets = appendNewFacets(zh, getLength(horizon_fcts_edges));
    initNewFacets(zh,point_index,new_facets,horizon_fcts,horizon_fcts_edges);
    //    printf("available points: ");
    //    printList(avail_points);
    //    printf("new facets : ");
    //    printList(new_facets);
    dividePointsBetweenNewFacets(zh, avail_points, new_facets, horizon_fcts);
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
    list_t available_points=emptyList();
    if (maxit>MAXIT)
        maxit=MAXIT;
    if (zh->pts.num_points!=0){
        zhullInitialFacets(zh);
        printZhull(zh);
        while((getLength(zh->facets_with_outsidepoints)>0)&&(cnt++<maxit))  {
            printf("//////////////// ITERATION %d ///////\n",cnt);
            fli%=getLength(zh->facets_with_outsidepoints);
            f=getFacetByIndex(zh->facets_with_outsidepoints,fli);
            pi=f->farthest_outside_point;
            removeVisibleFacetsGetHorizonAndAvailablePoints(zh,pi,f,
                    &horizon_fcts, &horizon_fcts_edges,
                    &available_points);
            removeValueFromList(&available_points, pi);
            makePyramidFacetsToHorizon(zh,pi,horizon_fcts,horizon_fcts_edges,
                    available_points);
            printZhull(zh);

            freeList(&horizon_fcts);
            freeList(&horizon_fcts_edges);
            freeList(&available_points);
            fli++;
        }
    }
}

void printZhull(const zhull_t * const zh) {
    index_t fi;
    list_t indices = emptyList();
    printf("zhull from %d points\n", zh->pts.num_points);
    printf("facets with outsidepoints: ");
    indices=findValueListInList(zh->facets_with_outsidepoints,zh->facets);
    printList(indices);
    freeList(&indices);
//    printf("facets with insidepoints: ");
//    printList(zh->facets_with_insidepoints);
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

