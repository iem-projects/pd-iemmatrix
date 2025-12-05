---
title: mtx_show
description: display a matrix
categories:
- object
pdcategory: General
see_also:
- mtx_print
- mtx_plot
inlets:
  1st:
  - type: matrix
    description: matrix to print to the console
outlets:
  1st:
  - type: matrix
    description: changed matrix
---

{{<pdmsg mtx_show >}} displays a matrix within the patch.
Elements can be edited via double click or by dragging.

Note, that this object is very slow and should only be used for small matrices.
