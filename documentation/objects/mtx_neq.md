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

$$C_{m\times n} = (A_{m\times n} \stackrel{?}{\not =}^\circ B_{m\times n}) \quad \equiv \quad c_{ij} = (a_{ij} \stackrel{?}{\not =}^\circ b_{ij})$$

$$C_{m\times n} = (A_{m\times n} \stackrel{?}{\not =}^\circ b) \quad \equiv \quad c_{ij} = (a_{ij} \stackrel{?}{\not =}^\circ b)$$

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
