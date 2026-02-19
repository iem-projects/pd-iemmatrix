---
title: mtx_matcher
description: match two sets of feature vectors
categories:
- object
pdcategory: Misc
see_also:
inlets:
  1st:
  - type: matrix
    description: set of vectors to match \(A_{n\times m}\)
  2nd:
  - type: matrix
    description: reference set of vectors \(B_{n\times m}\)
outlets:
  1st:
  - type: matrix
    description: permutated matrix \(A'_{n\times m}\)
  2nd:
  - type: matrix
    description: permutation matrix \(T_{n\times n}\)
---

Given two sets of feature vectors arranged as a matrix (each matrix row is a feature vector).

{{<pdobj mtx_matcher >}} will attempt to permutate the input set of vectors \(A\) in such a way,
that the euclidean distance between the permutated matrix \(A'\) and the reference set of vectors \(B\)
becomes minimal.

It can be used to for simple tracking algorithms,
where features vectors (e.g. coordinates)
are retrieved out-of-order (e.g. they are sorted by some arbitrary feature, like the x-Axis).
Matching the current set of feature vectors with a prediction of the feature vectors
will give sort the current set in such a way that it's row indices match the prediction.
(A very simply prediction algorithm for slowly-changing features is to just use the previous set.)

The permutation matrix can be used to transform the input set of vectors into the output set of vectors:

$$
A' = T \cdot A
$$


The result is ***not* guaranteed** to have an absolute minimum distance to the reference set of features.
