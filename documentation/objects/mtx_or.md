---
title: mtx_||
description: matrix logic OR
aliases:
- mtx_or
categories:
- object
pdcategory: Element Math
see_also:
- mtx_|
- mtx_&&
inlets:
  1st:
  - type: matrix
    description: input matrix A
  2nd:
  - type: matrix
    description: input matrix B
outlets:
  1st:
  - type: matrix
    description: boolean output matrix (A.||B)
---
Each element of the input matrices is converted to boolean value,
and afterwards the corresponding values are added.
The output is a boolean matrix (only ones and zeros).

$$C_{m,n} = A_{m,n} \lor ^{\circ} B_{m,n} \quad \equiv \quad c_{i,j} = a_{i,j} \lor  b_{i,j}$$
