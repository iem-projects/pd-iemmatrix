#ifndef QHULL_VECTOR_H
#define QHULL_VECTOR_H 
#include "list.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
 *  vector operations for zhull
 *
 * Copyright (c) 2012, Franz Zotter
 * IEM, Graz, Austria
 * 
 *
 */

typedef struct vec_ {
    float c[3];
} vector_t;

typedef struct points_ {
    vector_t *v;
    size_t num_points;
} points_t;

typedef struct plane_ {
    vector_t normal;
    vector_t point;
} plane_t;

typedef struct line_ {
    vector_t direction;
    vector_t point;
} line_t;


vector_t initVector (const float x, const float y, const float z);
float lengthVector(const vector_t v);
vector_t normalizeVector(const vector_t v);

plane_t initPlane (vector_t normal, const vector_t point);
line_t initLine (vector_t direction, const vector_t point);

points_t initPoints (const float *x,
        const float *y, const float *z, 
        const size_t num_points);
void freePoints (points_t *points);
size_t getNumPoints(const points_t points);
vector_t getPoint(const points_t points, const index_t index);

vector_t crossProduct (const vector_t v1, const vector_t v2);
float innerProduct (const vector_t v1, const vector_t v2);
float distancePointPoint(const vector_t a, const vector_t b);
float distancePointPlane (const vector_t point, const plane_t plane);
float distancePointLine (const vector_t point, const line_t line);
float distancePointLineOnPlane (const vector_t point, 
        const line_t line, const plane_t plane);
vector_t addVectors(const vector_t v1, const vector_t v2);
vector_t subtractVectors(const vector_t v1, const vector_t v2);
vector_t scaleVector(vector_t v1, const float f);
/*vector_t averagePoints(points_t points);
vector_t normalOfPoints(points_t points);
plane_t planeFromPoints(points_t points);*/
plane_t planeFromThreePoints (const vector_t a, const vector_t b, const vector_t c);
line_t lineFromTwoPoints (const vector_t a, const vector_t b);

vector_t averageListedPoints(const points_t points, const list_t list);
vector_t normalOfListedPoints(const points_t points, const list_t list);
vector_t directionOfListedPoints(const points_t points, const list_t list);
plane_t planeFromListedPoints(const points_t points, const list_t list);
line_t lineFromListedPoints(const points_t points, const list_t list);

void printPlane(const plane_t p);
void printLine(const line_t l);
void printVector(const vector_t v);
#endif /* QHULL_VECTOR_H */
