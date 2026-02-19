---
title: mtx_spherical_harmonics_rotator
description: create a rotation matrix for spherical harmonics
categories:
- object
pdcategory: Spherical Harmonics
see_also:
- mtx_spherical_harmonics
arguments:
- type: float
  description: order \(n\) of the spherical hamronics
inlets:
  1st:
  - type: list <float> <float> <float>
    description: Euler angles \(\alpha\), \(\beta\) & \(\gamma\) (in radian) for ZYZ rotation
outlets:
  1st:
  - type: matrix
    description: \(R_{{(n+1)^2}\times{(n+1)^2}}\) rotation matrix for spherical harmonics of order \(n\)
---

Creates a 3D rotation matrix for spherical harmonics of order \(n\).


Useful for manipulating Ambisonic encoded soundfields.
