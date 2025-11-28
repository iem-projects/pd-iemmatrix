---
title: mtx_eig
description: eigenvalue decomposition of a matrix
categories:
- object
pdcategory: Matrix Math
see_also:
- mtx_svd
inlets:
  1st:
  - type: matrix
    description: \(n\times n\) matrix \(\boldsymbol A\) to decompose
outlets:
  1st:
  - type: list
    description: list of real parts of the \(n\) eigenvalues \(\boldsymbol \sigma\)
  2nd:
  - type: list
    description: list of imaginary parts of the \(n\) eigenvalues \(\boldsymbol \sigma\)
  3:
  - type: matrix
    description: \(n\times n\) real part of the eigenvector matrix \(\boldsymbol V\)
  4:
  - type: matrix
    description: \(n\times n\) imaginary part of the eigenvector matrix \(\boldsymbol V\)
draft: false
---
This object depends on the GNU Scientific Library. The matrices lists involved have the properties:

$$ \boldsymbol A = \boldsymbol V \mathrm{diag}\\{\boldsymbol \sigma\\} \boldsymbol V^{-1} $$


