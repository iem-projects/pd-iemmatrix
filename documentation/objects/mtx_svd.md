---
title: mtx_svd
description: singular value decomposition of a matrix
categories:
- object
pdcategory: Matrix Math
see_also:
- mtx_qr
inlets:
  1st:
  - type: matrix
    description: \(n\times m\) matrix \(\boldsymbol A\) to decompose
outlets:
  1st:
  - type: matrix
    description: \(n\times n\) matrix of left-singular vectors \(\boldsymbol U\)
  1:
  - type: list
    description: list of singular values \(\boldsymbol s\)
  2:
  - type: matrix
    description: \(m\times m\) matrix of right-singular vectors \(\boldsymbol V\)
draft: false
---
This object depends on the GNU Scientific Library. The matrices lists involved have the properties:

$$ \boldsymbol A_{n\times m} = \boldsymbol U \mathrm{diag}_{n\times m}\\{\boldsymbol s\\} \boldsymbol{V}^\intercal $$

$$ \boldsymbol U\boldsymbol U^\intercal=\boldsymbol U^\intercal\boldsymbol U=\boldsymbol I_{n\times n} $$

$$ \boldsymbol V\boldsymbol V^\intercal=\boldsymbol V^\intercal\boldsymbol V=\boldsymbol I_{m\times m} $$

$$ \boldsymbol s_{\mathrm{min}\\{n,m\\}\times 1}=[s_1,s_2,\dots]^\intercal $$ with $$ s_1\geq s_2\geq s_3 \dots \geq 0 $$

