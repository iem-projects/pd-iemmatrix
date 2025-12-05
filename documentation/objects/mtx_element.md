---
title: mtx_element
description: set elements of a matrix
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_row
- mtx_col
inlets:
  1st:
  - type: float
    description: new value for selected matrix element(s)
  - type: matrix
    description: new values for all matrix elements
  2nd:
  - type: float float
    description: element 2d position
outlets:
  1st:
  - type: matrix
    description: modified matrix
arguments:
  - type: float float [float float]
    description: "initial matrix dimensions; optional: initial element position"
---

This object allows you to change specific elements of a matrix.
The element to change is selected by its 2D index to the 2nd outlet.

Indices are *`1`-based*, so the position {{<pdmsg 1 1 >}} refers to the element at the 1st column in the 1st row.

You can use `0` as a wildcard index, so the position {{<pdmsg 2 0 >}} refers to all elements in the 2nd row.
