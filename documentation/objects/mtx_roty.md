---
title: mtx_roty
description: 3d cartesian rotation matrix around Y-axis
categories:
- object
pdcategory: Coordinates
see_also:
- mtx_rot
- mtx_rotx
- mtx_rotz
- mtx_rotxyz
inlets:
  1st:
  - type: matrix
    description: rotation angle \(\theta\) (in radian)
outlets:
  1st:
  - type: matrix
    description: 3D rotation matrix
---

$$R_y = \begin{pmatrix}
\cos \theta & 0 & -\sin \theta \\
0 & 1 & 0 \\
\sin \theta & 0 & \cos \theta
\end{pmatrix}
$$
