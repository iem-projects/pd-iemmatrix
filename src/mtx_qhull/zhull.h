
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
    list_t facets;
    list_t facets_with_outsidepoints;
    list_t facets_with_insidepoints;
} zhull_t;

void calculateZHull(zhull_t *zh,int maxit);
void printZhull(const zhull_t * const  zh);
void freeZhull(zhull_t *zh);
zhull_t zhullInitPoints(const float *x, const float *y, 
        const float *z, const size_t num_points);
void printFacet(const zhull_t * const zh, 
      const facet_t * const f);
