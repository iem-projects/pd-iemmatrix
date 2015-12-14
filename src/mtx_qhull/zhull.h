#ifndef QHULL_ZHULL_H
#define QHULL_ZHULL_H

#include "vectors.h"
#include "list.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define TOL_OUTSIDEPOINT 1e-7
#define TOL_INSIDEPOINT  -1e-7
#define TOL_DEGENERATE  1e-6
#define MAXIT 1000000
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

typedef struct facet_ {
  plane_t plane;
  list_t corners;
  list_t outsideset;
  list_t insideset;
  size_t farthest_outside_point;
  list_t neighbors;
  float maxdistance;
} facet_t;

typedef struct zhull_ {
  points_t pts;
  list_t used_pts;
  list_t facets;
  list_t facets_with_outsidepoints;
  list_t facets_with_insidepoints;
} zhull_t;

int calculateZHull(zhull_t *zh);
index_t getTriangleCorner(const zhull_t * const zh,
                          const index_t triangle_idx,
                          const index_t corner_idx);
void printZhull(const zhull_t * const  zh);
void freeZhull(zhull_t *zh);
zhull_t zhullInitPoints(const float *x, const float *y,
                        const float *z, const size_t num_points);
void printFacet(const zhull_t * const zh,
                const facet_t * const f);
#endif /* QHULL_ZHULL_H */
