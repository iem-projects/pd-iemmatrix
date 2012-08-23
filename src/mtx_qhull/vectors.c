#include "vectors.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


vector_t initVector (float x, float y, float z) {
    vector_t vec={x, y, z};
    return vec;
}

vector_t normalizeVector(vector_t v) {
    float r=sqrtf(v.c[0]*v.c[0]+v.c[1]*v.c[1]+v.c[2]*v.c[2]);
    v.c[0]/=r;
    v.c[1]/=r;
    v.c[2]/=r;
    return v;
};

plane_t initPlane (vector_t normal, vector_t point) {
    plane_t plane;
    plane.point = point;
    plane.normal = normalizeVector(normal);
    return plane;
}

points_t initPoints (const float *x,
        const float *y, const float *z, 
        size_t num_points) {
    points_t points;
    size_t i;

    points.v = (vector_t*) malloc (sizeof(vector_t)*num_points);
    if (points.v!=0) {
//        printf("created %li\n",points.v);
        points.num_points = num_points;
        for (i=0; i<num_points; i++) {
            points.v[i] = initVector(x[i],y[i],z[i]);
        }
    }
    return points; 
}

void freePoints (points_t *points) {
//    printf("deleting %li\n",points->v);
    if (points!=0) {
        if (points->v!=0)
            free(points->v);
        points->v = 0;
        points->num_points = 0;
    }
}

vector_t crossProduct (vector_t v1, vector_t v2) {
    vector_t cp;
    cp.c[0]= v1.c[1]*v2.c[2]-v1.c[2]*v2.c[1];
    cp.c[1]=-v1.c[0]*v2.c[2]+v1.c[2]*v2.c[0];
    cp.c[2]= v1.c[0]*v2.c[1]-v1.c[1]*v2.c[0];
    return cp;
}

float innerProduct (vector_t v1, vector_t v2) {
    return v1.c[0]*v2.c[0] + v1.c[1]*v2.c[1] + v1.c[2]*v2.c[2];
}

float distancePointPlane (vector_t point, plane_t plane) {
    return innerProduct(point, plane.normal) - 
        innerProduct(plane.point, plane.normal);
}

vector_t addVectors(vector_t v1, vector_t v2) {
    vector_t v3;
    v3.c[0]=v1.c[0]+v2.c[0];
    v3.c[1]=v1.c[1]+v2.c[1];
    v3.c[2]=v1.c[2]+v2.c[2];
    return v3;
}

vector_t subtractVectors(vector_t v1, vector_t v2) {
    vector_t v3;
    v3.c[0]=v1.c[0]-v2.c[0];
    v3.c[1]=v1.c[1]-v2.c[1];
    v3.c[2]=v1.c[2]-v2.c[2];
    return v3;
}

vector_t scaleVector(vector_t v1, float f) {
    vector_t v2;
    v2.c[0]=f*v1.c[0];
    v2.c[1]=f*v1.c[1];
    v2.c[2]=f*v1.c[2];
    return v2;
}

vector_t averagePoints(points_t points) {
    vector_t m = initVector(0.0f, 0.0f, 0.0f);
    size_t i;
    for (i=0; i<points.num_points; i++)
        m=addVectors(points.v[i],m);
    m=scaleVector(m,1.0f/((float)points.num_points));
    return m;
}

vector_t normalOfPoints(points_t points) {
    vector_t m = averagePoints(points);
    vector_t n = initVector(0.0f, 0.0f, 0.0f);
    size_t i;
    for (i=0; i<points.num_points; i++) {
        n=addVectors(crossProduct(points.v[i],m),n);
    }
    return n;
}

plane_t planeFromPoints(points_t points) {
    vector_t p=averagePoints(points);
    vector_t n=normalOfPoints(points);
    plane_t pl=initPlane(n,p);
    return pl;
}

void printPlane(plane_t p) {
   printf("n=");
   printVector(p.normal);
   printf(", p=");
   printVector(p.point);
}

void printVector(vector_t v) {
   printf("[%5.2f,%5.2f,%5.2f], ", v.c[0],v.c[1],v.c[2]);
}

