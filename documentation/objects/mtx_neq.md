---
title: mtx_!=
description: element-wise comparison for non-equality of two matrices
categories:
- object
pdcategory: Element Math
aliases:
- mtx_neq
see_also:
- mtx_isequal
- "mtx_=="
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

$$C_{m,n} = (A_{m,n} \stackrel{?}{\not =}^\circ B_{m,n}) \quad \equiv \quad c_{i,j} = (a_{i,j} \stackrel{?}{\not =}^\circ b_{i,j})$$

$$C_{m,n} = (A_{m,n} \stackrel{?}{\not =}^\circ b) \quad \equiv \quad c_{i,j} = (a_{i,j} \stackrel{?}{\not =}^\circ b)$$

$$
[\begin{pmatrix}
1 & 0 \cr
0 & 2
\end{pmatrix} \stackrel{?}{\not =}^\circ \begin{pmatrix}
1 & 0 \cr
2 & 2
\end{pmatrix}] = \begin{pmatrix}
0 & 0 \cr
1 & 0
\end{pmatrix}
$$
