% A_fhp.m
%
% Takes measurements in the upper half-plane of the 2D Fourier transform.
%
% Usage: b = A_fhp(x, OMEGA)
%
% x - N vector
%
% b - K vector = [mean; real part(OMEGA); imag part(OMEGA)]
%
% OMEGA - K/2-1 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
%
% Written by: Justin Romberg, Caltech
% Created: October 2005
% Email: jrom@acm.caltech.edu
%

function y = A_fhp(x, OMEGA)

n = round(sqrt(length(x)));

yc = 1/n*fft2(reshape(x,n,n));
y = [yc(1,1); sqrt(2)*real(yc(OMEGA)); sqrt(2)*imag(yc(OMEGA))];

