---
title: mtx_convolver~
description: multiple-input-multiple-output (MIMO) convolvution matrix
categories:
- object
pdcategory: Audio Math
see_also:
inlets:
  1st:
  - type: array3
    description: contains the \(n\times m\times len\) impulse responses of filters connecting every of the \(m\) inputs to every of the \(n\)  outputs
  m:
  - type: signal
    description: multiple input signals that are to be convolved using uniformly partitioned overlap save with the FFT size \(2L\), where the partition length \(L\) that subdivides the total impulse response length \(len\) is controlled by the subpatch [block~] setting in Pd. The object permits amplitude-complementary or power-complementary time-varying updates within an L samples output crossfade
outlets:
  m:
  - type: signal
    description: multiple output signals of which each is obtained by filtering and summing with the m filters for every input of length L, for the respective output
arguments:
  - type: <symbol>
    description: \(-\mathrm{m}\) option to switch to multichannel input and output connections replacing the individual ins and outs
  - type: <symbol>
    description: \(\mathrm{pow}\) pow option to switch from amplitude-complementary to power-complementary crossfades during time-varying updates 
  - type: <symbol>
    description: file name to read array3 impulse response configuration from

draft: false
---

## array3
The array3 message is shaped as [array3 \<n\> \<m\> \<len\> \<ir-samples\> ( with the linear index 
$$s+i\times len+o\times m\times len, \qquad \text{indices: } s\dots\text{sample, } i\dots\text{input, } o\dots\text{output}$$ to shape the impulse response (ir)  sample sequence.

## time-varying MIMO convolver 
The input signal \\(x_i[s]\\) of the input index \\(i\\) and sample index \\(s\\) is divided by Pd into blocks \\(x_b[s]=x[bL+s]\\) of the length \\(L\\) with the local sample index \\(s=0\dots L\\).  Together with the previous block, it is being \\(2L\\) fast Fourier transformed for real-valued signals, at every DSP cycle indexed by \\(b\\):
$$
 X_{i,b}[k]=\mathrm{rFFT_{2L}}\\{[x_{b-1}[0]\dots x_{b-1}[L-1]], [x_{b}[0]\dots x_{b}[L-1]] \\}.
$$
The impulse responses \\(h_{i,o,c}[s]\\) for the inlet \\(i=0\dots m-1\\), outlet \\(o=0\dots n-1\\) and current/old \\(c=0,1\\) crossfade index are partitioned by `[mtx_convolver~]` into blocks of also \\(L\\), using the partition index \\(p=0\dots P-1\\) with \\(P=\lceil\mathrm{len}/L\rceil\\) and Fourier transformed after zero-padding to the size \\(2L\\), and the frequency bin index \\(k=0\dots L\\) is relevant for real-valued signals:
$$
 H_{i,o,p,c}[k]=\mathrm{rFFT_{2L}}\\{[h_{i,o,c}[pL],\dots h_{i,o,c}[pL+L-1]], [0,\dots,0]\\}.
$$
Cyclic convolution for the output \\(o\\) is obtained via complex-valued multiply accumulate for the relevant frequency bins \\(k=0\dots L\\) across the inputs\\(i=0\dots m-1\\) and partitions \\(p=0\dots P-1\\) for a current/old response set \\(c=0,1\\)
$$
  Y_{o,c}[k]=
             \sum_{i=0}^{m-1}
             \sum_{p=0}^{P-1} 
             X_{i,b-p}[k] H_{i,o,p,c}[k].
$$
For every output and current/old, the inverse fast Fourier transform for real-valued signals delivers
$$
  y_{o,c}[s]=\mathrm{riFFT_{2L}}\\{Y_{o,c}[k]\\},
$$
which, when there is an update, delivers for every output \\(o\\) the linear-convolution output of the length \\(L\\) from crossfading with an \\(L\\) point fade-out window \\(w[n]\\), which starts from the rIFFT sample \\(L\\)
$$
  y_{o}[s]=w[s]\thinspace y_{o,c}[s+L]+(1-w[s])\thinspace y_{o,(c+1)\\%2}[s+L],
$$
or with the option \\(\mathrm{pow}\\) for power-complementary output crossfade:
$$
  y_{o}[s]=\sqrt{w[s]}\thinspace y_{o,c}[s+L]+\sqrt{1-w[s]}\thinspace y_{o,(c+1)\\%2}[s+L],
$$
or if there is no update
$$
  y_{o}[s]=\thinspace y_{o,0}[s+L].
$$
The window \\(w[s]\\) is implemented as a Hann window of the length \\(L\\).
