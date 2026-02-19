---
title: mtx_pack~
description: pack signal vectors into a matrix
categories:
- object
pdcategory: Audio Math
see_also:
- mtx_unpack~
arguments:
  - type: float
    description: number of signal inlets
inlets:
  n:
  - type: signal
    description: input signals (possibly multichannel)
outlets:
  1st:
  - type: matrix
    description: signal vectors
  2nd:
  - type: channels <float>
    description: number of discrete input channels (number of rows in the output matrix)
  - type: blocksize <float>
    description: blocksize of the input signals (number of columns in the output matrix)
  - type: dimen <float> <float>
    description: number of channels and blocksize in a single message
---

Takes the current block of all discrete input signals and outputs them as a {{<pdmsg matrix >}} message
(making it possible to apply matrix math on the signals).

Converting the matrix back to a signal can be done with {{<pdobj "mtx_unpack~" >}}.

Note that converting a signal to a message and back again will introduce one block delay.
