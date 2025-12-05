---
title: mtx_plot
description: plot a matrix
categories:
- object
pdcategory: General
see_also:
- mtx_print
- mtx_show
inlets:
  1st:
  - type: matrix
    description: matrix to display
---

{{<pdmsg mtx_plot >}} displays a matrix as a 2D image.
Values are silently clamped to the \\([0, +1]\\) range.

The default greyscale display can be colorised via a mapping table.

Note, that this object is very slow and should only be used for small matrices.
