---
title: mtx_rotx
description: 3d cartesian rotation matrix around X-axis
categories:
- object
pdcategory: Coordinates
see_also:
- mtx_rot
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
    description: 3D rotation matrix
---

$$R_x = \begin{pmatrix}
1 & 0 & 0 \\
0 & \cos \theta & -\sin \theta \\
0 & \sin \theta & \cos \theta
\end{pmatrix}
$$
