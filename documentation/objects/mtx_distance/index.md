---
title: mtx_distance
description: calculate euclidean distance matrix between two sets of vectors
categories:
- object
pdcategory: Misc
see_also:
- mtx_distance2
inlets:
  1st:
  - type: matrix
    description: LxN matrix
  2nd:
  - type: matrix
    description: MxN matrix
outlets:
  1st:
  - type: matrix
    description: LxM matrix
---

{{< pdobj mtx_distance >}} gets all the euclidean distances between two sets of vectors.


The euclidean distance between two vectors (of equal dimensions \\(N\\)) \\(a\\) and \\(b\\) is:
$$
\| \vec a -\vec b \| = \sqrt{\textstyle\sum_{n=1}^N (a_n-b_n)^2} = \sqrt{(a_1-b_1)^2 + (a_2-b_2)^2 + \dots}
$$


Assuming \\(A_{L\times N}\\) consists of \\(L\\) vectors of dimension \\(N\\)
and \\(B_{M\times N}\\) consists of \\(M\\) vectors of dimension \\(N\\),
then {{< pdobj mtx_distance >}} returns a matrix \\(C_{L\times M}\\),
where each row vector \\(\vec c_l = C_{l,\*}\\) corresponds to the euclidean distances between
the same-index row vector \\(\vec a_l = A_{l,\*}\\) and all the row vectors in \\(B\\):

$$
C = \begin{pmatrix}
\| \vec a_1 - \vec b_1 \| & \| \vec a_1 - \vec b_2 \| & \| \vec a_1 - \vec b_3 \| & \dots \\
\| \vec a_2 - \vec b_1 \| & \| \vec a_2 - \vec b_2 \| & \| \vec a_2 - \vec b_3 \| & \dots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}
$$


## Example

![distance between two matrices](mtx_distance.pd.svg)

The first row of the output matrix shows the distances between the reference vector
\\(\begin{pmatrix}0 & 0 & 0\end{pmatrix}\\) and the three vectors
\\(\begin{pmatrix}1 & 0 & 0\end{pmatrix}\\),
\\(\begin{pmatrix}0 & 1 & 0\end{pmatrix}\\), and
\\(\begin{pmatrix}0 & 0 & 1\end{pmatrix}\\) resp.
The distance is \\(1\\) in all cases.

The second row of the output matrix uses \\(\begin{pmatrix}1 & 2 & 3\end{pmatrix}\\) as the reference vector,
leading to the results
\\(\sqrt{0^2+2^2+3^2} = \sqrt{13} = 3.605\\),
\\(\sqrt{1^2+1^2+3^2} = \sqrt{11} = 3.316\\), and
\\(\sqrt{1^2+2^2+2^2} = \sqrt{9} = 3\\) resp.

## Performance considerations
If you are mainly interested in sorting distances between points,
you might consider using {{< pdobj mtx_distance2 >}} instead,
which gives you the *squares* of each distance
(the \\(\sqrt{\dots}\\) operation is typically considered rather costly).
