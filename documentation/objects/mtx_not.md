---
title: mtx_!
description: matrix logic NOT
aliases:
- mtx_not
categories:
- object
pdcategory: Element Math
see_also:
- mtx_||
- mtx_&&
inlets:
  1st:
  - type: matrix
    description: input matrix
outlets:
  1st:
  - type: matrix
    description: boolean output matrix (.!A)
---
Each element of the input matrix is converted to a boolean value,
and afterwards inverted.
The output is a boolean matrix (only ones and zeros).



$$B_{m,n} = \lnot A_{m,n} \quad \equiv \quad b_{i,j} = \lnot a_{i,j}$$
