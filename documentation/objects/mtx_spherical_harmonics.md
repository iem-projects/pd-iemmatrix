---
title: mtx_spherical_harmonics
description: evaluate the spherical harmonics as matrix
categories:
- object
pdcategory: Audio Math
see_also:
- mtx_circular_harmonics
inlets:
  1st:
  - type: matrix
    description: \(2\times L\) matrix with azimuth angles (row 1) and  zenith angles (row 2) in radians for all \(L\) directions to evaluate
outlets:
  1st:
  - type: matrix
    description: \(L\times (N+1)^2\) matrix with the spherical harmonics evaluated from zeroth to \(N^\mathrm{th}\) order using the numerical creation argument \(N\) of the object
draft: false
---
Useful for Ambisonic encoding and decoding. The object permits several different normalizations N3D (default), N3D4PI, SN3D, SN3D
