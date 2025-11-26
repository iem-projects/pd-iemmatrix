---
title: mtx_.cmul
description: element-wise multiplication of complex-valued matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_.cabs2
- mtx_.cdiv
- mtx_cabs2
- mtx_cinverse
- mtx_cmul
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
    description: real part of element-wise multiplication of (A+jB).*(C+jD)
  2nd:
  - type: matrix
    description: imaginary part of element-wise multiplication of (A+jB).*(C+jD)
---
