---
title: mtx_unpack~
description: unpack a matrix into signal vectors
categories:
- object
pdcategory: Audio Math
see_also:
- mtx_pack~
arguments:
  - type: float
    description: number of signal inlets
  - type: <symbol>
    description: "`-m` option to enable multichannel mode"
inlets:
  1st:
  - type: matrix
    description: signal vectors
outlets:
  n:
  - type: signal
    description: output signal; in multichannel mode, there's only a single outlet
---

Takes a matrix with signal vectors with *rows* channels and a blocksize of *columns*
(as created by  {{<pdobj "mtx_pack~" >}}),
and turns it into proper signals.

In non-multichannel mode, the number of output channels is determined by the argument.
In multichannel mode (`-m`), there is only a single outlet for a multichannel signal. The number of channels
in the multichannel signal depends on number of rows in the input matrix.
Due to the way Pd works, the number of channels in a multichannel signal cannot change once the DSP is turned on.
If you need to change it, restart the DSP.
