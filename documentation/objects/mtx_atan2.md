---
title: mtx_atan2
description: arctangent (atan) of matrix elements (four-quadrant version)
categories:
- object
pdcategory: Element Math
see_also:
- mtx_sin
- mtx_cos
- mtx_tan
- mtx_atan
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
$$C_{m,n} = \operatorname{atan2^\circ} (A_{m,n}, B_{m,n}) \quad \equiv \quad c_{i,j} = \operatorname{atan}(a_{i,j} / b_{i,j})$$
