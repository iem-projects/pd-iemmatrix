#include "vectors.h"


vector_t initVector (float x, float y, float z) {
    vector_t vec={x, y, z};
    return vec;
}

float lengthVector(vector_t v) {
    return sqrtf(v.c[0]*v.c[0]+v.c[1]*v.c[1]+v.c[2]*v.c[2]);
}
vector_t normalizeVector(vector_t v) {
    float r=lengthVector(v);
    v.c[0]/=r;
    v.c[1]/=r;
    v.c[2]/=r;
    return v;
}

plane_t initPlane (vector_t normal, vector_t point) {
    plane_t plane;
    plane.point = point;
    plane.normal = normalizeVector(normal);
    return plane;
}

line_t initLine (vector_t direction, vector_t point) {
    line_t line;
    line.point = point;
    line.direction = normalizeVector(direction);
    return line;
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

float distancePointPoint(const vector_t a,
        const vector_t b) {
    vector_t d=subtractVectors(b,a);
    return lengthVector(d);
}

float distancePointPlane (vector_t point, plane_t plane) {
    return innerProduct(point, plane.normal) - 
        innerProduct(plane.point, plane.normal);
}

float distancePointLine (vector_t point, line_t line) {
    return lengthVector(crossProduct(line.direction,
                subtractVectors(point,line.point)));
}

float distancePointLineOnPlane (vector_t point, 
        line_t line, plane_t plane) {
    vector_t normal_in_plane = normalizeVector(crossProduct(
                line.direction, plane.normal));
    return innerProduct(subtractVectors(point,line.point),normal_in_plane);
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

/*
vector_t averagePoints(points_t points) {
    vector_t m = initVector(0.0f, 0.0f, 0.0f);
    size_t i;
    for (i=0; i<getNumPoints(points); i++)
        m=addVectors(getPoint(points,m));
    m=scaleVector(m,1.0f/((float)getNumPoints(points)));
    return m;
}
vector_t normalOfPoints(points_t points) {
    vector_t m = averagePoints(points);
    vector_t n = initVector(0.0f, 0.0f, 0.0f);
    size_t i;
    for (i=0; i<getNumPoints(points); i++) {
        n=addVectors(crossProduct(getPoint(points,i),m),n);
    }
    return n;
}
plane_t planeFromPoints(points_t points) {
    vector_t p=averagePoints(points);
    vector_t n=normalOfPoints(points);
    plane_t pl=initPlane(n,p);
    return pl;
}
*/

plane_t planeFromThreePoints(const vector_t a, 
        const vector_t b, const vector_t c) {
    vector_t ab = subtractVectors(b,a);
    vector_t ac = subtractVectors(c,a);
    vector_t n = normalizeVector(crossProduct(ab,ac));
    return initPlane(n,a);
}

line_t lineFromTwoPoints(const vector_t a, 
        const vector_t b) {
    vector_t direction = subtractVectors(b,a);
    return initLine(direction,a);
}

size_t getNumPoints(const points_t points) {
    return points.num_points;
}

vector_t getPoint(const points_t points, 
        const index_t index) {
    if ((index>=0)&&(index<getNumPoints(points)))
        return points.v[index];
    return initVector(0.0f,0.0f,0.0f);
}

vector_t averageListedPoints(const points_t points, const list_t list) {
    index_t i;
    vector_t m = initVector(0.0f, 0.0f, 0.0f);
    for (i=0; i<getLength(list); i++) {
        m=addVectors(getPoint(points,getEntry(list,i)),m);
    }
    m=scaleVector(m,1.0f/((float)getLength(list)));
    return m;
}

vector_t normalOfListedPoints(const points_t points, const list_t list) {
    vector_t m = averageListedPoints(points,list);
    vector_t n = initVector(0.0f, 0.0f, 0.0f);
    vector_t d1,d2,c;
    index_t i;
    for (i=1; i<=getLength(list); i++) {
        d1=subtractVectors(getPoint(points,getEntry(list,i-1)),m);
        d2=subtractVectors(getPoint(points,getEntry(list,i%getLength(list))),m);
        c=crossProduct(d1,d2);
        n=addVectors(c,n);
    }
    return n;
}

vector_t directionOfListedPoints(const points_t points, const list_t list) {
    vector_t m = averageListedPoints(points,list);
    vector_t dir = initVector(0.0f, 0.0f, 0.0f);
    vector_t d,c;
    index_t i;
    for (i=1; i<getLength(list); i++) {
        d=subtractVectors(getPoint(points,getEntry(list,i-1)),
                getPoint(points,getEntry(list,i)));
        dir=addVectors(d,dir);
    }
    return dir;
}

plane_t planeFromListedPoints(const points_t points, const list_t list) {
    vector_t p=averageListedPoints(points,list);
    vector_t n=normalOfListedPoints(points,list);
    plane_t pl=initPlane(n,p);
    return pl;
}



line_t lineFromListedPoints(const points_t points, const list_t list) {
    vector_t p=averageListedPoints(points,list);
    vector_t n=directionOfListedPoints(points,list);
    line_t l=initLine(n,p);
    return l;
}

void printVector(vector_t v) {
   printf("[%5.2f,%5.2f,%5.2f], ", v.c[0],v.c[1],v.c[2]);
}

void printPlane(plane_t p) {
   printf("n=");
   printVector(p.normal);
   printf(", p=");
   printVector(p.point);
}

void printLine(line_t l) {
   printf("d=");
   printVector(l.direction);
   printf(", p=");
   printVector(l.point);
}


