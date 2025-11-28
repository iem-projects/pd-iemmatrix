---
title: matrix
description: store/create/manipulate a 2d matrix
categories:
- object
pdcategory: General
aliases:
- mtx
see_also:
- mtx_ones
- mtx_zeros
- mtx_diag
- mtx_diegg
- mtx_eye
- mtx_egg
- mtx_col
- mtx_element
- mtx_row
inlets:
  1st:
  - type: matrix
    description: store and output the matrix
  - type: bang
    description: output the stored matrix
  2nd:
  - type: matrix
    description: store the matrix
outlets:
  1st:
  - type: matrix
    description: the stored matrix
# FIXXME: arguments <symbol> is not a second argument but an alternative
# we need a better way to express this
arguments:
  - type: <float> <float>
    description: number of rows & number of columns
  - type: <symbol>
    description: file name to read matrix from

---

## matrix store
{{< pdobj mtx >}} stores a single matrix, just like {{< pdobj float >}} stores a single number.

You can also save the stored matrix to a file via the message `[write <fname>(`
to the matrix, or load a matrix via `[read <fname>(` from a file.


## matrix creation
It can also be used to create special matrices using the following messages:
| message | result                   | related object           |
|---------|--------------------------|--------------------------|
| `zeros` | all zeros                | {{< pdobj mtx_zeros >}}  |
| `ones`  | all ones                 | {{< pdobj mtx_ones >}}   |
| `eye`   | identity matrix          | {{< pdobj mtx_eye >}}    |
| `egg`   | backward identity matrix | {{< pdobj mtx_egg >}}    |
| `diag`  | diagonal matrix          | {{< pdobj mtx_diag >}}   |
| `diegg` | backward diagonal matrix | {{< pdobj mtx_diegg >}}  |
| `size`  | silently resize matrix   | {{< pdobj mtx_resize >}} |


## matrix manipulation
Finally you can query and change the contents of a matrix:

| message   | result                  | related object            |
|-----------|-------------------------|---------------------------|
| `row`     | get/set *n*th row       | {{< pdobj mtx_row >}}     |
| `col`     | get/set *m*th column    | {{< pdobj mtx_col >}}     |
| `element` | get/set element *(n,m)* | {{< pdobj mtx_element >}} |


# the `[matrix(` message
The {{< pdobj matrix >}} object understands (and emits) `matrix` messages.

This message describes a 2d matrix and has the form

```
matrix <int:rows> <int:columns> <list[float]:entries>
```

The *entries* are a list of \\(rows * columns\\) float atoms in row-major order.

E.g.
$$
\texttt{[matrix\quad 3 2\quad 10 20 \thinspace 30 40 \thinspace 50 60(}
\to
\begin{pmatrix}
10 & 20 \cr
30 & 40 \cr
50 & 60
\end{pmatrix}
$$
