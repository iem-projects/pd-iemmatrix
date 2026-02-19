---
title: mtx_rot
description: 2d cartesian rotation matrix
categories:
- object
pdcategory: Coordinates
see_also:
- mtx_rotx
- mtx_roty
- mtx_rotz
- mtx_rotxyz
inlets:
  1st:
  - type: float
    description: rotation angle \(\theta\) (in radian)
outlets:
  1st:
  - type: matrix
    description: 2D rotation matrix
---

$$R = \begin{pmatrix}
\cos \theta & -\sin \theta \\
\sin \theta & \cos \theta
\end{pmatrix}
$$
