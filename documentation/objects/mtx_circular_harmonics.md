---
title: mtx_circular_harmonics
description: evaluate the circular harmonics as matrix
categories:
- object
pdcategory: Spherical Harmonics
see_also:
- mtx_spherical_harmonics
inlets:
  1st:
  - type: matrix
    description: \(1\times L\) matrix with azimuth angles in radians for all \(L\) directions to evaluate
outlets:
  1st:
  - type: matrix
    description: \(L\times (2N+1)\) matrix with the circular harmonics evaluated from zeroth to \(N^\mathrm{th}\) order using the numerical creation argument \(N\) of the object
draft: false
---
Useful for 2D Ambisonic encoding and decoding. The object permits several different normalizations N2D (default), N2D2PI, SN2D
