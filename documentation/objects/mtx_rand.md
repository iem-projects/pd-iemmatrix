---
title: mtx_rand
description: create a matrix with random values
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_ones
inlets:
  1st:
  - type: matrix
    description: create an equally-sized matrix
  - type: float
    description: create a square matrix with the given dimensions
  - type: float float
    description: create a matrix with the given dimensions
  - type: seed
    description: set the seed for the pseudo random generator
outlets:
  1st:
  - type: matrix
    description: matrix with random elements
---

Generated pseudo random values are in the range of \\(\left[ 0, 1 \right) \\).
