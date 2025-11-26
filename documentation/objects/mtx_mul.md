---
title: mtx_*
description: matrix multiplication
categories:
- object
pdcategory: Matrix Math
aliases:
- mtx_mul
see_also:
- matrix
- mtx_+
- mtx_-
- mtx_.*
- mtx_./
- mtx_.^
inlets:
  1st:
  - type: matrix
    description: left-hand operand
  2nd:
  - type: matrix
    description: right-hand operand
  - type: float
    description: right-hand operand
outlets:
  1st:
  - type: matrix
    description: result
---


$$C_{m,p} = A_{m,n} \cdot B_{n,p}$$

$$C_{m,n} = A_{m,n} \cdot b \quad \equiv \quad c_{i,j} = a_{i,j} \cdot b $$
