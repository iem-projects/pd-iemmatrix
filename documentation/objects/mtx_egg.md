---
title: mtx_egg
description: create an exchange matrix
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_egg
- mtx_diegg
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
    description: backward identity matrix
---


$$
I_{m,n} = \begin{pmatrix}
0 & \dots & 0 & 1\cr
0 & \dots & 1 & 0 \cr
\vdots & \vdots & \ddots & \vdots \cr
1 & 0 & \dots & 0
\end{pmatrix}
$$
