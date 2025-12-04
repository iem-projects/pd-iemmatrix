---
title: mtx_roll
description: shift columns of a matrix
categories:
- object
pdcategory: Matrix Transformation
see_also:
- mtx_scroll
inlets:
  1st:
  - type: matrix
    description: input matrix
  2nd:
  - type: float
    description: shift amount
outlets:
  1st:
  - type: matrix
    description: output matrix
---

Shifts matrix columns by *shift amount* to the right, wrapping around.
If *shift amount* is negative, columns are shifted to the left.

$$
\operatorname{roll}(\begin{pmatrix}
1 & 0 & 0 & 0 & 0 \cr
0 & 2 & 0 & 0 & 0 \cr
0 & 0 & 3 & 0 & 0 \cr
0 & 0 & 0 & 4 & 0 \cr
0 & 0 & 0 & 0 & 5
\end{pmatrix}, 2) = \begin{pmatrix}
0 & 0 & 1 & 0 & 0 \cr
0 & 0 & 0 & 2 & 0 \cr
0 & 0 & 0 & 0 & 3 \cr
4 & 0 & 0 & 0 & 0 \cr
0 & 5 & 0 & 0 & 0
\end{pmatrix}
$$
