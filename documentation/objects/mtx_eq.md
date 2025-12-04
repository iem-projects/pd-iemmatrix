---
title: mtx_==
description: element-wise comparison for equality of two matrices
categories:
- object
pdcategory: Element Math
aliases:
- mtx_eq
see_also:
- mtx_isequal
- "mtx_!="
- "mtx_>"
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

$$C_{m,n} = (A_{m,n} \stackrel{?}{=}^\circ B_{m,n}) \quad \equiv \quad c_{i,j} = (a_{i,j} \stackrel{?}{=} b_{i,j})$$

$$C_{m,n} = (A_{m,n} \stackrel{?}{=}^\circ b) \quad \equiv \quad c_{i,j} = (a_{i,j} \stackrel{?}{=} b)$$

## Examples

$$
[\begin{pmatrix}
1 & 0 \cr
0 & 2
\end{pmatrix} \stackrel{?}{=}^\circ \begin{pmatrix}
1 & 0 \cr
2 & 2
\end{pmatrix}] = \begin{pmatrix}
1 & 1 \cr
0 & 1
\end{pmatrix}
$$
