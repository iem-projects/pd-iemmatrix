---
title: mtx_inverse
description: inverse (or pseudo-inverse) of a matrix
categories:
- object
pdcategory: Matrix Math
see_also:
inlets:
  1st:
  - type: matrix
    description: \(n\times n\) input matrix
outlets:
  1st:
  - type: matrix
    description: \(n\times n\) output matrix
---

Calculates the inverse of a matrix:
$$ Y_{n\times n} = X_{n\times n}^{-1} $$

so that \\(Y_{n\times n} *  X_{n\times n} = I\\).

If the input matrix  \\(X_{n\times m}\\) is not square (that is: \\(n != m\\)),
{{<pdobj mtx_inverse >}} will automatically calculate the *pseudo-inverse* instead:


$$
\begin{cases}
Y_{m\times n} = (X^\mathsf{T} \cdot X)^{-1} \cdot X^\mathsf{T} ,  & \text{if $n > m$} \cr
Y_{m\times n} = X^\mathsf{T} \cdot (X \cdot X^\mathsf{T})^{-1} ,  & \text{if $n < m$}
\end{cases}
$$
