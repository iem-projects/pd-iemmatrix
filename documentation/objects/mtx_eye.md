---
title: mtx_eye
description: create an identity matrix
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_egg
- mtx_diag
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
    description: identity matrix
---


$$
I_{m,n} = \begin{pmatrix}
1 & 0 & \dots & 0 \cr
0 & 1 & \dots & 0 \cr
\vdots & \vdots & \ddots & \vdots \cr
0  & 0 & \dots  & 1
\end{pmatrix}
$$
