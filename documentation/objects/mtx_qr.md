---
title: mtx_qr
description: calculate QR decomposition of a matrix
categories:
- object
pdcategory: Matrix Transformation
see_also:
- mtx_svd
inlets:
  1st:
  - type: matrix
    description: \(n\times m\) matrix \(\boldsymbol A\) to decompose
outlets:
  1st:
  - type: matrix
    description: \(n\times n\) matrix \(\boldsymbol Q\)
  2nd:
  - type: matrix
    description: \(n\times m\) matrix \(\boldsymbol R\)
draft: false
---
This object depends on the GNU Scientific Library. The matrices lists involved have the properties:

$$ \boldsymbol A_{n\times m} = \boldsymbol Q \boldsymbol R $$

$$ \boldsymbol Q^\intercal\boldsymbol Q=\boldsymbol Q\boldsymbol Q^\intercal =\bm I $$

$$ \boldsymbol R = \mathrm{triu}_{n\times m}\{\dots\} $$
