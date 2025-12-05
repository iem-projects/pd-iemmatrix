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


$$C_{m\times p} = A_{m\times n} \cdot B_{n\times p}$$

$$C_{m\times n} = A_{m\times n} \cdot b \quad \equiv \quad c_{ij} = a_{ij} \cdot b $$
