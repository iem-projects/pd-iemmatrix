---
title: mtx_./
description: element-wise matrix division
categories:
- object
pdcategory: Element Math
aliases:
- mtx_div
see_also:
- matrix
- mtx_+
- mtx_-
- mtx_.*
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



$$C_{m,n} = A_{m,n} \oslash B_{m,n} \quad \equiv \quad c_{i,j} = a_{i,j} / b_{i,j} $$

$$C_{m,n} = A_{m,n} \oslash b \quad \equiv \quad c_{i,j} = a_{i,j} / b $$
