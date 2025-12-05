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
$$C_{m\times n} = (A_{m\times n}+jB_{m\times n})^{|\circ|} \quad \equiv \quad c_{hi} = |a_{hi}+j{b_{hi}}| = \sqrt{a^2_{hi}+b^2_{hi}}$$
