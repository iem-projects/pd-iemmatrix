---
title: mtx_max2
description: get element-wise maxima of two matrices
categories:
- object
pdcategory: Element Math
see_also:
- mtx_min2
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
$$C_{m,n} = {\max}^\circ (A_{m,n}, B_{m,n}) \quad \equiv \quad c_{i,j} = \max(a_{i,j}, b_{i,j})$$
