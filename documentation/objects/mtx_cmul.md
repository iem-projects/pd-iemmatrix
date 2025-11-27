---
title: mtx_cmul
description: multiplication of complex-valued matrix
categories:
- object
pdcategory: Matrix Math
see_also:
- mtx_*
- mtx_.cmul
inlets:
  1st:
  - type: matrix
    description: real part of 1st input matrix (A)
  2nd:
  - type: matrix
    description: imaginary part of 1st input matrix (B)
  3rd:
  - type: matrix
    description: real part of 2nd input matrix (C)
  4th:
  - type: matrix
    description: imaginary part of 2nd input matrix (D)
outlets:
  1st:
  - type: matrix
    description: X - real part of multiplication of (A+jB)(C+jD)
  2nd:
  - type: matrix
    description: Y - imaginary part of multiplication of (A+jB)(C+jD)
---

$$
Z_{m,n} = (A_{m,n}+jB_{m,n}) \cdot  (C_{m,n}+jD_{m,n}) \\\\
X_{m,n} = \Re ({Z_{m,n}}) \\\\
Y_{m,n} = \Im ({Z_{m,n}})
$$
