---
title: mtx_isequal
description: compare a matrix for equality with another (scalar result)
categories:
- object
pdcategory: Misc
see_also:
- "mtx_=="
- "mtx_!="
- "mtx_>"
- "mtx_>="
- "mtx_<"
inlets:
  1st:
  - type: matrix
    description: left-hand operand \(A_{m\times n}\)
  2nd:
  - type: matrix
    description: right-hand operand \(B_{m_n}\)
  - type: float
    description: right-hand operand \(b\)
outlets:
  1st:
  - type: matrix
    description: result
---

Returns `1` if all elements of the input matrix \\(A_{m\times n}\\) are equal to the corresponding elements of \\(B_{m\times n}\\) (resp. \\(b\\)), and otherwise `0`.


$$
c = (A_{m\times n} \stackrel{?}{=} B_{m\times n}) =
\begin{cases}
1 & \text{if } a_{ij} = b_{ij} \forall i,j, \cr
0 & \text{otherwise.}
\end{cases}
$$


$$
c = (A_{m\times n} \stackrel{?}{=} b) =
\begin{cases}
1 & \text{if } a_{ij} = b \forall i,j, \cr
0 & \text{otherwise.}
\end{cases}
$$

## Examples

$$
[\begin{pmatrix}
1 & 0 \cr
0 & 2
\end{pmatrix} \stackrel{?}{=} \begin{pmatrix}
1 & 0 \cr
0 & 2
\end{pmatrix}] = 1
$$

$$
[\begin{pmatrix}
1 & 0 \cr
0 & 2
\end{pmatrix} \stackrel{?}{=} \begin{pmatrix}
1 & 1 \cr
2 & 2
\end{pmatrix}] = 0
$$

