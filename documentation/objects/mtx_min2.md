---
title: mtx_min2
description: get element-wise minima of two matrices
categories:
- object
pdcategory: Element Math
see_also:
- mtx_max2
- mtx_minmax
inlets:
  1st:
  - type: matrix
    description: right hand operand
  2nd:
  - type: matrix
    description: right hand operand
outlets:
  1st:
  - type: matrix
    description: output
---
$$C_{m,n} = {\min}^\circ (A_{m,n}, B_{m,n}) \quad \equiv \quad c_{i,j} = \min(a_{i,j}, b_{i,j})$$
