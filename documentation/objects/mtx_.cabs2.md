---
title: mtx_.cabs2
description: element-wise absolute value of complex-valued matrix
categories:
- object
pdcategory: Element Math
see_also:
- mtx_.cdiv
- mtx_.cmul
- mtx_cabs2
- mtx_cinverse
- mtx_cmul
inlets:
  1st:
  - type: matrix
    description: real part of input matrix (A)
  2nd:
  - type: matrix
    description: imaginary part of input matrix (B)
outlets:
  1st:
  - type: matrix
    description: absolute values of A+jB
---
$$C_{m,n} = (A_{m,n}+jB_{m,n})^{|\circ|} \quad \equiv \quad c_{h,i} = |a_{h,i}+j{b_{h,i}}| = \sqrt{a^2_{h,i}+b^2_{h,i}}$$
