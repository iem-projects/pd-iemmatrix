---
title: matrix~
description: matrix multiplication of signals (with interpolation)
categories:
- object
pdcategory: Obsolete
see_also:
- mtx_*~
inlets:
  0..M-1:
  - type: signal
    description: input signals
  M:
  - type: matrix
    description: matrix to multiply input signals with
  M+1:
  - type: float
    description: interpolation time
outlets:
  0..N-1:
  - type: signal
    description: output signals
---

This object is obsolete.
It shouldn't be used in new code anymore.
Please use {{< pdobj "mtx_*~" >}} instead.
