---
title: mtx_dbtopow
description: calculate power value from dB for each matrix element
categories:
- object
pdcategory: Element Math
see_also:
- mtx_dbtorms
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

Applies {{< pdobj dbtopow >}} to each matrix element.
