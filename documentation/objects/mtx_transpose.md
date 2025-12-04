---
title: mtx_transpose
description: transpose a matrix
categories:
- object
pdcategory: Matrix Transformation
see_also:
inlets:
  1st:
  - type: matrix
    description: \(n\times m\) input matrix
outlets:
  1st:
  - type: matrix
    description: \(m\times n\) output matrix
---

$$ Y_{m\times n} = X_{n\times m}^\mathsf{T} $$
