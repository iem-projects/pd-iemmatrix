---
title: mtx_cinverse
description: matrix inverse of complex-valued matrix
categories:
- object
pdcategory: Matrix Math
see_also:
- mtx_inverse
inlets:
  1st:
  - type: matrix
    description: real part of \(n\times n\) input matrix
  2nd:
  - type: matrix
    description: imaginary part of \(n\times n\) input matrix
outlets:
  1st:
  - type: matrix
    description: real part of \(n\times n\) output matrix
  2nd:
  - type: matrix
    description: imaginary part of \(n\times n\) output matrix
---

Calculates the inverse of a complex-valued matrix:
$$ Y_{n\times n} = X_{n\times n}^{-1} $$
