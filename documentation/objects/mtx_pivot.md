---
title: mtx_pivot
description: pivot transform a matrix
categories:
- object
pdcategory: Matrix Math
see_also:
inlets:
  1st:
  - type: matrix
    description: \(n\times m\) input matrix
outlets:
  1st:
  - type: matrix
    description: \(n\times m\) output matrix
  2nd:
  - type: matrix
    description: \(n\times n\) pre-multiplication matrix
  3rd:
  - type: matrix
    description: \(m\times m\) post-multiplication matrix
---

This will transform the columns and rows of the input matrix \\(X\\), so that the result \\(Y\\) will have "all" maximum values in the diagonal.
The maximum of the matrix will be located at the upper-left corner.

The first outlet is the pivot-transformed matrix.
The other outlets are the elementary matrices \\(A\\) and \\(B\\)
that have to be pre-multiplied (row-transform) resp. post-multiplied (column-tranform)
with the original matrix to get the pivot-tranformation.

$$
Y_{n\times m} = A_{n\times n} \cdot X_{n\times m} \cdot B_{m\times m}
$$

This is useful for things like de-pivoting.


## Example

$$
X_{m,n} = \begin{pmatrix}
1 & 2 & 3 \cr
4 & 5 & 6 \cr
7 & 8 & 9 \cr
10 & 11 & 12
\end{pmatrix}
$$

$$
Y_{m,n} = \begin{pmatrix}
12 & 11 & 10 \cr
9 & 8 & 7 \cr
6 & 5 & 4 \cr
3 & 2 & 1 
\end{pmatrix}
= \begin{pmatrix}
0 & 0 & 0 & 1 \cr
0 & 0 & 1 & 0 \cr
0 & 1 & 0 & 0 \cr
1 & 0 & 0 & 0
\end{pmatrix} \cdot \begin{pmatrix}
1 & 2 & 3 \cr
4 & 5 & 6 \cr
7 & 8 & 9 \cr
10 & 11 & 12
\end{pmatrix} \cdot \begin{pmatrix}
0 & 0 & 1 \cr
0 & 1 & 0 \cr
1 & 0 & 0
\end{pmatrix}
$$
