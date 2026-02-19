---
title: mtx_cart2sph
description: convert cartesian to spherical coordinates (matrix version)
categories:
- object
pdcategory: Coordinates
see_also:
- mtx_sph2cart
inlets:
  1st:
  - type: matrix
    description: \(3\times L\) matrix containing x (row 1), y (row 2), z (row 3) of the \(L\) cartesian coordinates to convert
outlets:
  1st:
  - type: matrix
    description: \(3\times L\) matrix containing radius (row 1), azimuth (row 2), zenith (row 3) for the resulting \(L\) coordinates
draft: false
---
