# mtx convolver



## mtx convolver 

This repo contains source code for real-time multiple-input-multiple-output convolution with time-variant changes (cross fades). The implementation uses uniformly partitioned overlap-save block-convolution and keeps the number of input and output ffts as small as possible (in+out, and not in x out, and only a common one for all partitions). The implementation does not destroy input buffers on update, but uses an initialized number of ins, outs, partitions, and two output buffers per output for temporal cross-fade. It is kept simplistic by not considering multi-processing and multi-threading, and to permit easy modification for research purposes.

A Pd example implementation demonstrates the real-time mimo convolver.