---
title: mtx_scroll
description: shift rows of a matrix
categories:
- object
pdcategory: Matrix Transformation
see_also:
- mtx_roll
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

Shifts matrix rows by *shift amount* to the bottom, wrapping around.
If *shift amount* is negative, rows are shifted to the top.

$$
\operatorname{scroll}(\begin{pmatrix}
1 & 0 & 0 & 0 & 0 \cr
0 & 2 & 0 & 0 & 0 \cr
0 & 0 & 3 & 0 & 0 \cr
0 & 0 & 0 & 4 & 0 \cr
0 & 0 & 0 & 0 & 5
\end{pmatrix}, 2) = \begin{pmatrix}
0 & 0 & 0 & 4 & 0 \cr
0 & 0 & 0 & 0 & 5 \cr
1 & 0 & 0 & 0 & 0 \cr
0 & 2 & 0 & 0 & 0 \cr
0 & 0 & 3 & 0 & 0
\end{pmatrix}
$$
