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



$$C_{m\times n} = A_{m\times n} \oslash B_{m\times n} \quad \equiv \quad c_{ij} = a_{ij} / b_{ij} $$

$$C_{m\times n} = A_{m\times n} \oslash b \quad \equiv \quad c_{ij} = a_{ij} / b $$
