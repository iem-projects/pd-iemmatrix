---
title: mtx_rifft
description: real-valued inverse FFT for each matrix row
categories:
- object
pdcategory: Misc
see_also:
- [mtx_ifft]
- [mtx_rfft]
- [mtx_fft]
inlets:
  1st:
  - type: matrix
    description: input matrix with nfft/2+1 columns for the bins \(k=0\dots\text{nfft}/2\)
outlets:
  1st:
  - type: matrix
    description: output matrix with nfft columns, and nfft=\(2^l\) for \(l=2,3,\dots\)

draft: false
---
