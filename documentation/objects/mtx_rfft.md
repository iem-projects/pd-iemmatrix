---
title: mtx_rfft
description: real-valued FFT for each matrix row
categories:
- object
pdcategory: Digital Fourier Transformation
see_also:
- mtx_fft
- mtx_ifft
- mtx_rifft
inlets:
  1st:
  - type: matrix
    description: input matrix with nfft=\(2^l\) columns for \(l=2,3,\dots\)
outlets:
  1st:
  - type: matrix
    description: output matrix with nfft/2+1 columns for the bins \(k=0\dots\text{nfft}/2\)
draft: false
---
