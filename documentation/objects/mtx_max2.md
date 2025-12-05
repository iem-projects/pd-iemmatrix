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
$$C_{m\times n} = {\max}^\circ (A_{m\times n}, B_{m\times n}) \quad \equiv \quad c_{ij} = \max(a_{ij}, b_{ij})$$
