---
title: mtx_int
description: compute integer value of each matrix element
categories:
- object
pdcategory: Element Math
see_also:
inlets:
  1st:
  - type: matrix
    description: input matrix
outlets:
  1st:
  - type: matrix
    description: output matrix
---
Truncates any decimals of the matrix elements, so we only have integers.

$$B_{m\times n} = \operatorname{int^\circ} (A_{m\times n}) \quad \equiv \quad b_{ij} = \operatorname{int}(a_{ij})$$
