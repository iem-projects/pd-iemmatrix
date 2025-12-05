---
title: mtx_row
description: set rows of a matrix
categories:
- object
pdcategory: Matrix Creation
see_also:
- matrix
- mtx_element
- mtx_col
inlets:
  1st:
  - type: list
    description: new values for selected matrix row(s)
  - type: matrix
    description: new values for all matrix elements of an \(A_{n\times m}\) matrix
  2nd:
  - type: float
    description: row index
outlets:
  1st:
  - type: matrix
    description: modified matrix
arguments:
  - type: float float [float]
    description: "initial matrix dimensions; optional: initial row index"
---

This object allows you to change a specific row of a matrix.
The row to change is selected by its index to the 2nd outlet.

Indices are *`1`-based*, so the index {{<pdmsg 2 >}} refers to the 2nd row.

You can use `0` as a wildcard index, so the index {{<pdmsg 0 >}} refers to all rows in the matrix.

If the length of the list of new values exceeds the matrix width, the additional values are discarded.
If the length of this list is exactly `1`, all elements of the row are set to this value.
Otherwise, the length of this list must not be \\(\lt m\\) (less than the matrix width).
