% A_f.m
%
% Takes "scrambled Fourier" measurements.
%
% Usage: b = A_f(x, OMEGA, P)
%
% x - N vector
%
% b - K vector = [real part; imag part]
%
% OMEGA - K/2 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
%
% P - Permutation to apply to the input vector.  Fourier coeffs of
%     x(P) are calculated.
%     Default = 1:N (no scrambling).
%
% Written by: Justin Romberg, Caltech
% Created: October 2005
% Email: jrom@acm.caltech.edu
%

function b = A_f(x, OMEGA, P)

N = length(x);
if (nargin < 3), P = 1:N;  end

fx = 1/sqrt(N)*fft(x(P));
b = [sqrt(2)*real(fx(OMEGA)); sqrt(2)*imag(fx(OMEGA))];

