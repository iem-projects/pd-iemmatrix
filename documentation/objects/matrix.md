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
  - type: "message [matrix <rows> <columns> <float entries>( defining the matrix 
           in the sequence i = column + (row-1) x columns"
    description:
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

You can also save the stored matrix to a file via the message [write <fname>( to the matrix, or load a matrix via [read <fname>( from a file.


## matrix creation
It can also be used to create special matrices using the following messages:
| message | result                   | related object          |
|---------|--------------------------|-------------------------|
| `zeros` | all zeros                | {{< pdobj mtx_zeros >}} |
| `ones`  | all ones                 | {{< pdobj mtx_ones >}}  |
| `eye`   | identity matrix          | {{< pdobj mtx_eye >}}   |
| `egg`   | backward identity matrix | {{< pdobj mtx_egg >}}   |
| `diag`  | diagonal matrix          | {{< pdobj mtx_diag >}}  |
| `diegg` | backward diagonal matrix | {{< pdobj mtx_diegg >}} |
| `size`  | silently resize matrix   |                         |


## matrix manipulation
Finally you can query and change the contents of a matrix:

| message   | result                  | related object            |
|-----------|-------------------------|---------------------------|
| `row`     | get/set *n*th row       | {{< pdobj mtx_row >}}     |
| `col`     | get/set *m*th column    | {{< pdobj mtx_col >}}     |
| `element` | get/set element *(n,m)* | {{< pdobj mtx_element >}} |
