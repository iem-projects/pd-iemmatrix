---
title: mtx_&&
description: matrix logic AND
aliases:
- mtx_and
categories:
- object
pdcategory: Element Math
see_also:
- mtx_||
- mtx_&
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
    description: boolean output matrix (A.&&B)
---
Each element of the input matrices is converted to boolean value,
and afterwards the corresponding values are multiplied.
The output is a boolean matrix (only ones and zeros).

$$C_{m,n} = A_{m,n} \land ^{\circ} B_{m,n} \quad \equiv \quad c_{i,j} = a_{i,j} \land  b_{i,j}$$
