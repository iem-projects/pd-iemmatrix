---
title: mtx_log
description: compute the natural logarithm of each element of a matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_exp
- mtx_.^
inlets:
  1st:
  - type: matrix
    description: input matrix
outlets:
  1st:
  - type: matrix
    description: output matrix
---

$$B_{m,n} = \ln^\circ (A_{m,n}) \quad \equiv \quad b_{i,j} = \ln(a_{i,j})$$
