---
title: mtx_rotz
description: 3d cartesian rotation matrix around Z-axis
categories:
- object
pdcategory: Coordinates
see_also:
- mtx_rot
- mtx_rotx
- mtx_roty
- mtx_rotxyz
inlets:
  1st:
  - type: float
    description: rotation angle \(\theta\) (in radian)
outlets:
  1st:
  - type: matrix
    description: 3D rotation matrix
---

$$R_z = \begin{pmatrix}
\cos \theta & \sin \theta & 0 \\
-\sin \theta & \cos \theta & 0 \\
0 & 0 & 1
\end{pmatrix}
$$
