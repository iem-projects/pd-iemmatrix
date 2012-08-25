#ifndef QHULL_VECTOR_H
#define QHULL_VECTOR_H 
#include "list.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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


vector_t initVector (float x, float y, float z);
float lengthVector(vector_t v);
vector_t normalizeVector(vector_t v);

plane_t initPlane (vector_t normal, vector_t point);
line_t initLine (vector_t direction, vector_t point);

points_t initPoints (const float *x,
        const float *y, const float *z, 
        size_t num_points);
void freePoints (points_t *points);
size_t getNumPoints(const points_t points);
vector_t getPoint(const points_t points, const index_t index);

vector_t crossProduct (vector_t v1, vector_t v2);
float innerProduct (vector_t v1, vector_t v2);
float distancePointPoint(const vector_t a, const vector_t b);
float distancePointPlane (vector_t point, plane_t plane);
float distancePointLine (vector_t point, line_t line);
float distancePointLineOnPlane (vector_t point, 
        line_t line, plane_t plane);
vector_t addVectors(vector_t v1, vector_t v2);
vector_t subtractVectors(vector_t v1, vector_t v2);
vector_t scaleVector(vector_t v1, float f);
/*vector_t averagePoints(points_t points);
vector_t normalOfPoints(points_t points);
plane_t planeFromPoints(points_t points);*/
plane_t planeFromThreePoints (vector_t a, vector_t b, vector_t c);
line_t lineFromTwoPoints (vector_t a, vector_t b);

vector_t averageListedPoints(points_t points, list_t list);
vector_t normalOfListedPoints(points_t points, list_t list);
vector_t directionOfListedPoints(points_t points, list_t list);
plane_t planeFromListedPoints(points_t points, list_t list);
line_t lineFromListedPoints(points_t points, list_t list);

void printPlane(plane_t p);
void printLine(line_t l);
void printVector(vector_t v);
#endif /* QHULL_VECTOR_H */
