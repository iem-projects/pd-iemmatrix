---
title: mtx_zeros
description: create a matrix with all elements set to 0
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_ones
inlets:
  1st:
  - type: matrix
    description: create an equally-sized matrix
  - type: float
    description: create a square matrix with the given dimensions
  - type: float float
    description: create a matrix with the given dimensions
outlets:
  1st:
  - type: matrix
    description: matrix with all elements 0
---

$$
0_{m,n} = \begin{pmatrix}
0 & 0 & \dots & 0 \cr
0 & 0 & \dots & 0 \cr
\vdots & \vdots & \ddots & \vdots \cr
0  & 0 & \dots  & 0
\end{pmatrix}
$$
