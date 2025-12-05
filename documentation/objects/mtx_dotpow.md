---
title: mtx_.^
description: calculate element-wise raising to the power of x
categories:
- object
pdcategory: Element Math
aliases:
- mtx_pow
see_also:
- matrix
- mtx_+
- mtx_-
- mtx_.*
- mtx_./
- mtx_*
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


$$C_{m\times n} = A_{m\times n}^{\circ B_{m\times n}} \quad \equiv \quad c_{ij} = a_{ij}^{b_{ij}} $$

$$C_{m\times n} = A_{m\times n}^{\circ b} \quad \equiv \quad c_{ij} = a_{ij}^b $$
