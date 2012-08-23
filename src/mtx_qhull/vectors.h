#ifndef QHULL_VECTORS_H
#define QHULL_VECTORS_H

#include <sys/types.h>
#include "list.h"

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

vector_t initVector (float x, float y, float z);
vector_t normalizeVector(vector_t v);

plane_t initPlane (vector_t normal, vector_t point);

points_t initPoints (const float *x,
        const float *y, const float *z, 
        size_t num_points);
void freePoints (points_t *points);

vector_t crossProduct (vector_t v1, vector_t v2);
float innerProduct (vector_t v1, vector_t v2);
float distancePointPlane (vector_t point, plane_t plane);
vector_t addVectors(vector_t v1, vector_t v2);
vector_t subtractVectors(vector_t v1, vector_t v2);
vector_t scaleVector(vector_t v1, float f);
vector_t averagePoints(points_t points);
vector_t normalOfPoints(points_t points);
plane_t planeFromPoints(points_t points);

vector_t averageListedPoints(points_t points, list_t list);
vector_t normalOfListedPoints(points_t points, list_t list);
plane_t planeFromListedPoints(points_t points, list_t list);

void printPlane(plane_t p);
void printVector(vector_t v);

#endif /* QHULL_VECTORS_H */
