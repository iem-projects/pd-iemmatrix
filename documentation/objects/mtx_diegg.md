---
title: mtx_diegg
description: extract anti diagonal from matrix & create anti diagonal matrix
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
    description: extract anti diagonal from matrix
  - type: list
    description: create a square anti diagonal matrix with elements as provided in list
outlets:
  1st:
  - type: list
    description: anti diagonal from input matrix as list
  - type: matrix
    description: anti diagonal matrix created from input list
---

$$
\begin{pmatrix}
 1 &  2 &  3 &  4 \cr
 5 &  6 &  7 &  8 \cr
 9 & 10 & 11 & 12 \cr
13 & 14 & 15 & 16
\end{pmatrix}
\to
\begin{bmatrix}
4\cr 7\cr 10\cr 13\cr
\end{bmatrix}
$$

$$
\begin{bmatrix}
1\cr 2\cr 4\cr 3\cr
\end{bmatrix}
\to
\begin{pmatrix}
0 & 0 & 0 & 1 \cr
0 & 0 & 2 & 0 \cr
0 & 4 & 0 & 0 \cr
3 & 0 & 0 & 0
\end{pmatrix}
$$
