---
title: mtx_ones
description: create a matrix with all elements set to 1
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_zeros
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
    description: matrix with all elements 1
---



$$
J_{m\times n} = \begin{pmatrix}
1 & 1 & \dots & 1 \cr
1 & 1 & \dots & 1 \cr
\vdots & \vdots & \ddots & \vdots \cr
1  & 1 & \dots  & 1
\end{pmatrix}
$$
