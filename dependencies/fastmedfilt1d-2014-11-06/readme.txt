This ZIP file contains a fast Matlab implementation of 1D median filtering.
With the MEX core routine compiled using a decent compiler, compared against
Matlab's own proprietary toolbox implementation, this algorithm easily
achieves 10:1 performance gains for large window sizes. Note that, although
there are more efficient algorithms used in 2D image processing, these are
restricted to integer-valued data.

If you use this code for your research, please cite [1].

References:

[1] M.A. Little, N.S. Jones (2010), Sparse Bayesian Step-Filtering for High-
Throughput Analysis of Molecular Machine Dynamics in 2010 IEEE International
Conference on Acoustics, Speech and Signal Processing, 2010. ICASSP 2010
Proceedings.: Dallas, TX, USA (in press)

ZIP file contents:

fastmedfilt1d.m - The main routine. This calls the MEX core routine described
 below. Ensure this file is placed in the same directory as the MEX files
 below. Typing 'help fastmedfilt1d' gives usage instructions.

fastmedfilt1d_core.c - Core routine for performing running median filtering,
 written in C with Matlab MEX integration.

fastmedfilt1d_core.mexw32 - Above code compiled as a Matlab version 7 or
 greater library for direct use with Matlab under Windows 32. Place the
 library in a directory accessible to Matlab and invoke as with any other
 function.
