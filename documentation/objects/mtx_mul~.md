---
title: mtx_*~
aliases:
- mtx_mul~
description: multiply signals with a matrix
categories:
- object
pdcategory: Audio Math
see_also:
inlets:
  1st:
  - type: matrix
    description: vector or matrix to multiply input signals wit
  m:
  - type: signal
    description: multiple input signals that are understood as column vector per time sampling instant
  last: 
  - type: float
    description: interpolation time in ms to ramp from old matrix entries to new ones without crackling
outlets:
  n:
  - type: signal
    description: output signals that are understood as output colunm vector per sample
arguments:
  - type: <float> <float> <float>
    description: number of outputs and number of inputs and interpolation time in ms
  - type: <symbol>
    description: -m option to switch to multichannel input and output connections replacing the individual ins and outs
draft: false
---
$$
\begin{pmatrix} y_1[i]\cr
\vdots\cr
y_n[i]
\end{pmatrix}=
\boldsymbol{A}
\begin{pmatrix} x_1[i]\cr
\vdots\cr
x_m[i]
\end{pmatrix}
$$
