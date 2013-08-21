% At_f.m
%
% Adjoint for "scrambled Fourier" measurements.
%
% Usage: x = At_f(b, N, OMEGA, P)
%
% b - K vector = [real part; imag part]
%
% N - length of output x
%
% OMEGA - K/2 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
%
% P - Permutation to apply to the input vector.  Fourier coeffs of
%     x(P) are embedded.
%     Default = 1:N (no scrambling).
%
% Written by: Justin Romberg, Caltech
% Created: October 2005
% Email: jrom@acm.caltech.edu
%


function x = At_f(b, N, OMEGA, P)

if (nargin < 4),  P = 1:N;  end

K = length(b);
fx = zeros(N,1);
fx(OMEGA) = sqrt(2)*b(1:K/2) + i*sqrt(2)*b(K/2+1:K);
x = zeros(N,1);
x(P) = sqrt(N)*real(ifft(fx));
