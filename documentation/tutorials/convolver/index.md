---
title: "MIMO convolution"
# the lower the weight, the earlier it is sorted
weight: 10
draft: false
---

{{<pdobj "mtx_convolver~">}} provides functionality for MIMO (*M*ultiple *I*nputs *M*ultiple *O*utputs) signal processing,
using uniformly partitioned block convolution for efficient FIR filtering.

## Basics

The {{<pdobj "mtx_*~">}} object simply multiplies a matrix \\(A\\) with the the current multi-channel input \\(\vec{x}\\) to obtain the multi-channel output \\(\vec{y}\\);
$$
\vec{y} = A \cdot \vec{x}
$$
which is
$$
y_i[n]=\sum_{j=1}^J a_{ij}x_j[n]
$$

Rather than applying a single weight \\(a_{ij}\\) to each input channel,
the {{<pdobj "mtx_convolver~">}} applies a convolution at each node:

$$
y_i[n] = \sum_{j=1}^J \mathit{IR}_{ij} * x_j[n]
$$

### the `array3` message
The complete set of (equal length) impulse responses (IRs) can be seen as a three-dimensional matrix \\(\mathit{IR}_{ijk}\\).

Since iemmatrix's `[matrix ...(` message is defined as a two-dimensional matrix, a new message `[array3 ...(` is used.

The format is
```
[array3 <outs> <ins> <len> <data...>(
```

with `outs` and `ins` being the number of output channels (matrix rows) resp. input channels (matrix columns),
and `len` being the length of each impulse-response (in samples).
The following `data` is the actual impulse responses as a list of floats (in row-major order).

E.g.
```
[array3  3 2 4  1 0 0 0  2 2 2 2  3 0 0 3  1.1 1.2 1.3 1.4  2.1 2.2 2.3 2.4 3.1 -3.2 3.3 -3.4(
```

Is a \\(3*2\\) matrix with 6 impulse responses (4 samples each):

| IR  | values    |
|-----|-----------|
| \\(\mathit{IR}_{1,1}\\) | \\((1, 0, 0, 0)\\) |
| \\(\mathit{IR}_{1,2}\\) | \\((2, 2, 2, 2)\\) |
| \\(\mathit{IR}_{1,3}\\) | \\((3, 0, 0, 3)\\) |
| \\(\mathit{IR}_{2,1}\\) | \\((1.1, 1.2, 1.3, 1.4)\\) |
| \\(\mathit{IR}_{2,2}\\) | \\((2.1, 2.2, 2.3, 2.4)\\) |
| \\(\mathit{IR}_{2,3}\\) | \\((3.1,-3.2, 3.3,-3.4)\\) |



### partition size

The partition size for the partitioned convolution can be controlled via
Pd's block size, using the built-in {{<pdobj "block~">}} object.

When changing the impulse response matrix,
{{<pdobj "mtx_convolver~">}} will interpolate between the new and the old IR
within a single DSP block (so the interpolation time is `blocksize` (in samples)).


## special basic cases
With special convolution matrices we can use {{<pdobj "mtx_convolver~">}} to
mimic the behaviour of {{<pdobj "mtx_*~">}}  resp.  {{<pdobj "partconv~">}}.


### SISO convolution
The traditional object in Pd for SISO (*S*ingle *I*nput *S*ingle *O*utput) convolution is {{<pdobj "partconv~">}} from the [bsaylor](https://deken.puredata.info/library/bsaylor) library:

![SISO convolution with [partconv~]](simple_partconv~.svg)


Using a  \\(1 * 1 * J \\) convolution matrix, we can achieve the same result with {{<pdobj "mtx_convolver~">}}.

![SISO convolution with [mtx_convolver~]](simple_mtx_convolver~.svg)

The main differences are:
- we pass the impulse response via a `matrix` rather than a table
- the partition size is set via {{<pdobj "block~">}} rather than an argument


### matrix multiplication

Using IRs of length `1`, we can mimic the {{<pdobj "mtx_*~">}} (with a fixed interpolation time):

![matrix multiplication](simple_multiplication.svg)

Note that because the multiplication is done as convolution,
{{<pdobj "mtx_convolver~">}} is much slower than {{<pdobj "mtx_*~">}}.


# 2x2 cross-feed reverb

![2x2 cross-feed reverb: overall structure](example_2x2_crossfeed.svg)

This patch implements a dynamic 2×2 cross-feed reverb
where the impulse responses are generated and updated in real-time during playback.


A synthesizer subpatch ({{< pdobj "pd blip-blop" >}}) generates a short stereo test tone.
This signal is fed into a subpatch containing {{< pdobj "mtx_convolver~" >}} with a the 2×2 convolution matrix.
The subpatch is configured with a large block size ({{< pdobj "block~ 4096" >}}) (for efficiency reasons).


The matrix uses four impulse responses: two short IRs for the direct paths and two longer, modulated IRs for the cross-feed paths.
These IRs are synthesized within the {{< pdobj "pd cross-reverb" >}} subpatch and sent to the convolver as a list message.
When the user modifies the reverb parameters, the IRs are regenerated and swapped.

![2x2 cross-feed reverb: dynamically create IRs](example_2x2_crossfeed_makeIRs.svg)

The left-hand side of the patch creates two unit impulses (so the left and right input channel will pass through unaltered on their resp. channel).
For the cross-over, two more IRs are generated (on the right-hand side),
by using a random burst with an exponential decay curve, based on the input parameter (middle)
that is modulated with a cosine wave (far right). The modulation of the right resp. left IRs use different period lengths.
