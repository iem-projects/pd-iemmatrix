---
title: mtx_atan
description: arctangent (atan) of matrix elements
categories:
- object
pdcategory: Element Math
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
$$B_{m\times n} = \operatorname{atan^\circ} (A_{m\times n}) \quad \equiv \quad b_{ij} = \operatorname{atan}(a_{ij})$$
