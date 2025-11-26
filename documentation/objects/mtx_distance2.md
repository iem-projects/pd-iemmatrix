---
title: mtx_distance2
description: calculate the matrix of squared euclidean distances between a set of vectors
categories:
- object
pdcategory: Misc
see_also:
- mtx_distance
inlets:
  1st:
  - type: matrix
    description: LxN matrix
  2nd:
  - type: matrix
    description: MxN matrix
outlets:
  1st:
  - type: matrix
    description: LxM matrix
---

{{< pdobj mtx_distance2 >}} gets all the *squared* euclidean distances between two sets of vectors.

It behaves identical to the {{< pdobj mtx_distance >}} object,
but does not perform the final \\(\sqrt{\dots}\\) operation, so it is a bit more efficient.
