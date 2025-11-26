---
title: mtx_powtodb
description: convert a power value to dB for each element of a matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_dbtopow
- mtx_dbtorms
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

Applies {{< pdobj powtodb >}} to each matrix element.
