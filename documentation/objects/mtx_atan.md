---
title: mtx_atan
description: arctangent (atan) of matrix elements
categories:
- object
pdcategory: Trigonometric Functions
see_also:
- mtx_sin
- mtx_cos
- mtx_tan
- mtx_atan2
inlets:
  1st:
  - type: matrix
    description: input
outlets:
  1st:
  - type: matrix
    description: output
---
$$B_{m,n} = \operatorname{atan^\circ} (A_{m,n}) \quad \equiv \quad b_{i,j} = \operatorname{atan}(a_{i,j})$$
