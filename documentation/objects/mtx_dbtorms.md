---
title: mtx_dbtorms
description: calculate linear value from dB for each matrix element
categories:
- object
pdcategory: Element Math
see_also:
- mtx_dbtopow
- mtx_powtodb
- mtx_rmstodb
inlets:
  1st:
  - type: matrix
    description: NxM input matrix
outlets:
  1st:
  - type: matrix
    description: NxM output matrix
---

Applies {{< pdobj dbtorms >}} to each matrix element.
