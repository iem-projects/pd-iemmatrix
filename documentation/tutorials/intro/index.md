---
title: "A quick introduction to iemmatrix"
# the lower the weight, the earlier it is sorted
weight: 1
draft: false
---

## hello world

A *matrix* is a two-dimensional array of numbers.

*iemmatrix* uses a special {{<pdmsg matrix >}} message to send matrices from one object to another.
This message consists of the selector `matrix` followed by the number of *rows* and *columns*, and finally the data in row-major order:

$$
\texttt{[matrix   2 3  10 20 30  40 50 60(}
\to
\begin{pmatrix}
10 & 20 & 30 \cr
40 & 50 & 60
\end{pmatrix}
$$


![first patch](1.pd.svg)

The actual message format usually does not need to concern the user,
as matrices are created and processed by specialized objects.



## creating special matrices

The {{< pdobj mtx >}} object understands messages to create some special matrices:


![create special matrices](2.pd.svg)

There's also a number of specialized objects to create these matrices, e.g.
{{< pdobj mtx_ones >}}, {{< pdobj mtx_zeros >}},
{{< pdobj mtx_eye >}}, {{< pdobj mtx_egg >}},
{{< pdobj mtx_diag >}}, {{< pdobj mtx_diegg >}}
and {{< pdobj mtx_rand >}}.


## binary operators
![binops](3.pd.svg)

The {{< pdobj "mtx_+" >}}, {{< pdobj "mtx_-" >}},
{{< pdobj "mtx_.*" >}}, {{< pdobj "mtx_./" >}},
{{< pdobj "mtx_.^" >}},
{{< pdobj "mtx_>" >}}, {{< pdobj "mtx_>=" >}},
{{< pdobj "mtx_<" >}}, {{< pdobj "mtx_<=" >}},
{{< pdobj "mtx_==" >}} & {{< pdobj "mtx_!=" >}}
{{< pdobj "mtx_&&" >}} and {{< pdobj "mtx_||" >}},
{{< pdobj "mtx_&" >}} and {{< pdobj "mtx_|" >}}
objects all perform element-wise operations.
The operations are the same as their Pd-vanilla non-matrix counterparts.


For proper matrix multiplication use {{< pdobj "mtx_*" >}}.


## basic manipulation

![basic manipulation](4.pd.svg)

## advanced manipulation

![advanced manipulation](5.pd.svg)


## Examples

### finding the maximum, replacing it with -777

![finding the maximum](6.pd.svg)


### lin-swept sine

![lin-swept sine](7.pd.svg)


### simple peak picker

![simple peak picker](8.pd.svg)
