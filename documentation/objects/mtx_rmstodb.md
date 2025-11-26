---
title: mtx_rmstodb
description: convert a linear value into dB for each element of a matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_dbtopow
- mtx_dbtorms
- mtx_powtodb
inlets:
  1st:
  - type: matrix
    description: NxM input matrix
outlets:
  1st:
  - type: matrix
    description: NxM output matrix
---

Applies {{< pdobj rmstodb >}} to each matrix element.
