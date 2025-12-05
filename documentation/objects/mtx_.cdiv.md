---
title: mtx_.cdiv
description: element-wise division of complex-valued matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_.cabs2
- mtx_.cmul
- mtx_cabs2
- mtx_cinverse
- mtx_cmul
- mtx_./
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
    description: real part of element-wise division of (A+jB)./(C+jD)
  2nd:
  - type: matrix
    description: imaginary part of element-wise division of (A+jB)./(C+jD)
---

$$X_{m\times n} = (A_{m\times n}+jB_{m\times n}) \oslash  (C_{m\times n}+jD_{m\times n}) \quad \equiv \quad x_{hi} = (a_{hi}+j{b_{hi}})/(c_{hi}+j{d_{hi}})$$
