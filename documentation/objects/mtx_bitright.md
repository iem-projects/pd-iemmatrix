---
title: mtx_>>
description: matrix logic bitwise rightshift
alias:
- mtx_bitright
categories:
- object
pdcategory: Element Math
see_also:
- mtx_>>
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
    description: boolean output matrix (A.>>B)
---
Each element of the input matrices is converted to an integer value.
Afterwards right signed bit shift is performed on the elements of matrix A.
The shift amount is read from the corresponding elements in matrix B.
