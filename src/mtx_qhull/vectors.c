#include "vectors.h"
/*
 *  vector operations for zhull
 *
 * Copyright (c) 2012, Franz Zotter
 * IEM, Graz, Austria
 *
 *
 */



vector_t initVector (const float x, const float y, const float z)
{
  vector_t vec= {x, y, z};
  return vec;
}

float lengthVector(const vector_t v)
{
  return sqrtf(v.c[0]*v.c[0]+v.c[1]*v.c[1]+v.c[2]*v.c[2]);
}
vector_t normalizeVector(vector_t v)
{
  float r=lengthVector(v);
  v.c[0]/=r;
  v.c[1]/=r;
  v.c[2]/=r;
  return v;
}

plane_t initPlane (vector_t normal, const vector_t point)
{
  plane_t plane;
  plane.point = point;
  plane.normal = normalizeVector(normal);
  return plane;
}

line_t initLine (vector_t direction, const vector_t point)
{
  line_t line;
  line.point = point;
  line.direction = normalizeVector(direction);
  return line;
}

points_t allocatePoints (const size_t num_points)
{
  points_t points;
  points.v = (vector_t *) malloc(sizeof(vector_t)*num_points);
  if (points.v!=0) {
    points.num_points = num_points;
  } else {
    points.num_points = 0;
  }
  return points;
}

points_t initPoints (const float *x,
                     const float *y, const float *z,
                     const size_t num_points)
{
  points_t points = allocatePoints(num_points);
  size_t i;

  for (i=0; i<getNumPoints(points); i++) {
    points.v[i] = initVector(x[i],y[i],z[i]);
  }
  return points;
}

void freePoints (points_t *points)
{
//    printf("deleting %li\n",points->v);
  if (points!=0) {
    if (points->v!=0) {
      free(points->v);
    }
    points->v = 0;
    points->num_points = 0;
  }
}

void reallocatePoints (points_t *points, const size_t num_points)
{
  if ((num_points>0)&&(points!=0)) {
    if (getNumPoints(*points)==0) {
      *points=allocatePoints(num_points);
    } else {
      points->v = (vector_t *) realloc(points->v,sizeof(vector_t)*num_points);
      if (points->v!=0) {
        points->num_points=num_points;
      } else {
        points->num_points=0;
      }
    }
    if (points->v!=0) {
      points->num_points = num_points;
    }
  } else {
    freePoints(points);
  }
}

void appendPoints(points_t *points,
                  const float *x, const float *y,
                  const float *z, const size_t num_points)
{
  const size_t n=getNumPoints(*points);
  size_t i,j;
  reallocatePoints(points,getNumPoints(*points)+num_points);
  for (i=n,j=0; i<getNumPoints(*points); i++, j++) {
    points->v[i] = initVector(x[j],y[j],z[j]);
  }
}

vector_t crossProduct (const vector_t v1, const vector_t v2)
{
  vector_t cp;
  cp.c[0]= v1.c[1]*v2.c[2]-v1.c[2]*v2.c[1];
  cp.c[1]=-v1.c[0]*v2.c[2]+v1.c[2]*v2.c[0];
  cp.c[2]= v1.c[0]*v2.c[1]-v1.c[1]*v2.c[0];
  return cp;
}

float innerProduct (const vector_t v1, const vector_t v2)
{
  return v1.c[0]*v2.c[0] + v1.c[1]*v2.c[1] + v1.c[2]*v2.c[2];
}

float distancePointPoint(const vector_t a,
                         const vector_t b)
{
  vector_t d=subtractVectors(b,a);
  return lengthVector(d);
}

float distancePointPlane (const vector_t point, const plane_t plane)
{
  return innerProduct(point, plane.normal) -
         innerProduct(plane.point, plane.normal);
}

float distancePointLine (const vector_t point, const line_t line)
{
  return lengthVector(crossProduct(line.direction,
                                   subtractVectors(point,line.point)));
}

float distancePointLineOnPlane (vector_t const point,
                                const line_t line, const plane_t plane)
{
  vector_t normal_in_plane = normalizeVector(crossProduct(
                               line.direction, plane.normal));
  return innerProduct(subtractVectors(point,line.point),normal_in_plane);
}

vector_t addVectors(const vector_t v1, const vector_t v2)
{
  vector_t v3;
  v3.c[0]=v1.c[0]+v2.c[0];
  v3.c[1]=v1.c[1]+v2.c[1];
  v3.c[2]=v1.c[2]+v2.c[2];
  return v3;
}

vector_t subtractVectors(const vector_t v1, const vector_t v2)
{
  vector_t v3;
  v3.c[0]=v1.c[0]-v2.c[0];
  v3.c[1]=v1.c[1]-v2.c[1];
  v3.c[2]=v1.c[2]-v2.c[2];
  return v3;
}

vector_t scaleVector(vector_t v1, const float f)
{
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
                             const vector_t b, const vector_t c)
{
  vector_t ab = subtractVectors(b,a);
  vector_t ac = subtractVectors(c,a);
  vector_t n = normalizeVector(crossProduct(ab,ac));
  return initPlane(n,a);
}

line_t lineFromTwoPoints(const vector_t a,
                         const vector_t b)
{
  vector_t direction = subtractVectors(b,a);
  return initLine(direction,a);
}

size_t getNumPoints(const points_t points)
{
  return points.num_points;
}

vector_t getPoint(const points_t points,
                  const index_t index)
{
  if ((index>=0)&&(index<getNumPoints(points))) {
    return points.v[index];
  }
  return initVector(0.0f,0.0f,0.0f);
}

vector_t averageListedPoints(const points_t points, const list_t list)
{
  index_t i;
  vector_t m = initVector(0.0f, 0.0f, 0.0f);
  for (i=0; i<getLength(list); i++) {
    entry_t e=getEntry(list,i);
    m=addVectors(getPoint(points,entry_getIndex(&e)),m);
  }
  m=scaleVector(m,1.0f/((float)getLength(list)));
  return m;
}

vector_t normalOfListedPoints(const points_t points, const list_t list)
{
  vector_t m = averageListedPoints(points,list);
  vector_t n = initVector(0.0f, 0.0f, 0.0f);
  vector_t d1,d2,c;
  index_t i;
  for (i=1; i<=getLength(list); i++) {
    entry_t e=getEntry(list,i-1);
    d1=subtractVectors(getPoint(points,entry_getIndex(&e)),m);
    e=getEntry(list,i%getLength(list));
    d2=subtractVectors(getPoint(points,entry_getIndex(&e)),m);
    c=crossProduct(d1,d2);
    n=addVectors(c,n);
  }
  return n;
}

vector_t directionOfListedPoints(const points_t points, const list_t list)
{
  vector_t m = averageListedPoints(points,list);
  vector_t dir = initVector(0.0f, 0.0f, 0.0f);
  vector_t d,c;
  index_t i;
  for (i=1; i<getLength(list); i++) {
    entry_t e1=getEntry(list,i-1);
    entry_t e2=getEntry(list,i  );

    d=subtractVectors(getPoint(points,entry_getIndex(&e1)),
                      getPoint(points,entry_getIndex(&e2)));
    dir=addVectors(d,dir);
  }
  return dir;
}

plane_t planeFromListedPoints(const points_t points, const list_t list)
{
  vector_t p=averageListedPoints(points,list);
  vector_t n=normalOfListedPoints(points,list);
  plane_t pl=initPlane(n,p);
  return pl;
}



line_t lineFromListedPoints(const points_t points, const list_t list)
{
  vector_t p=averageListedPoints(points,list);
  vector_t n=directionOfListedPoints(points,list);
  line_t l=initLine(n,p);
  return l;
}

void printVector(const vector_t v)
{
  printf("[%5.2f,%5.2f,%5.2f], ", v.c[0],v.c[1],v.c[2]);
}

void printPlane(const plane_t p)
{
  printf("n=");
  printVector(p.normal);
  printf(", p=");
  printVector(p.point);
}

void printLine(const line_t l)
{
  printf("d=");
  printVector(l.direction);
  printf(", p=");
  printVector(l.point);
}
