---
title: mtx_|
description: matrix bitwise OR
aliases:
- mtx_bitor
categories:
- object
pdcategory: Element Math
see_also:
- mtx_||
- mtx_&
inlets:
  1st:
  - type: matrix
    description: input matrix A
  2nd:
  - type: matrix
    description: input matrix B
outlets:
  1st:
  - type: matrix
    description: output matrix (A.|B)
---
Each element of the input matrices is converted to integer and interpreted as a bitfield.
Afterwards an "or" is performed on the corresponding bits within the
corresponding matrix elements.

$$C_{m\times n} = A_{m\times n} | ^{\circ} B_{m\times n} \quad \equiv \quad c_{ij} = a_{ij} | b_{ij}$$
