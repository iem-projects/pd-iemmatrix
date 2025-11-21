---
title: mtx_diag
description: extract diagonal from matrix & create diagonal matrix
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_eye
- mtx_diegg
inlets:
  1st:
  - type: matrix
    description: extract diagonal from matrix
  - type: list
    description: create a square diagonal matrix with elements as provided in list
outlets:
  1st:
  - type: list
    description: diagonal from input matrix as list
  - type: matrix
    description: diagonal matrix created from input list
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
1\cr 6\cr 11\cr 16\cr
\end{bmatrix}
$$

$$
\begin{bmatrix}
1\cr 2\cr 4\cr 3\cr
\end{bmatrix}
\to
\begin{pmatrix}
1 & 0 & 0 & 0 \cr
0 & 2 & 0 & 0 \cr
0 & 0 & 4 & 0 \cr
0 & 0 & 0 & 3
\end{pmatrix}
$$
