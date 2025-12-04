---
title: mtx_&
description: matrix bitwise AND
aliases:
- mtx_bitand
categories:
- object
pdcategory: Element Math
see_also:
- mtx_&&
- mtx_|
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
    description: output matrix (A.&B)
---
Each element of the input matrices is converted to integer and interpreted as a bitfield.
Afterwards an "and" is performed on the corresponding bits within the
corresponding matrix elements.


$$C_{m,n} = A_{m,n} \And ^{\circ} B_{m,n} \quad \equiv \quad c_{i,j} = a_{i,j} \And  b_{i,j}$$
