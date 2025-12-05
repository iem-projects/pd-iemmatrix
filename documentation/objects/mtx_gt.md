---
title: mtx_>
description: element-wise comparison for strict superiority of two matrices
categories:
- object
pdcategory: Element Math
aliases:
- mtx_gt
see_also:
- mtx_isequal
- "mtx_=="
- "mtx_!="
- "mtx_>="
- "mtx_<"
- "mtx_<="
inlets:
  1st:
  - type: matrix
    description: left-hand operand
  2nd:
  - type: matrix
    description: right-hand operand
  - type: float
    description: right-hand operand
outlets:
  1st:
  - type: matrix
    description: result
---

Returns `1` for each element of \\(A\\) that is greater than the corresponding element in \\(B\\).

$$C_{m\times n} = (A_{m\times n} \stackrel{?}{\gt}^\circ B_{m\times n}) \quad \equiv \quad c_{ij} = (a_{ij} \stackrel{?}{\gt} b_{ij})$$

$$C_{m\times n} = (A_{m\times n} \stackrel{?}{\gt}^\circ b) \quad \equiv \quad c_{ij} = (a_{ij} \stackrel{?}{\gt} b)$$

## Examples

$$
[\begin{pmatrix}
1 & 1 \cr
0 & 2
\end{pmatrix} \stackrel{?}{\gt}^\circ \begin{pmatrix}
1 & 0 \cr
2 & 2
\end{pmatrix}] = \begin{pmatrix}
0 & 1 \cr
0 & 0
\end{pmatrix}
$$
