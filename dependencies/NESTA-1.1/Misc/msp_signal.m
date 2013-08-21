%
% [x0,b,A,At,supind,Omega,b_exact]=msp_signal(N,M,K,U,Ut,Dyna,NoiseLevel,Mode,Repeat)
%
%   It generates K-sparse signals or loads wavelet coefficients from a file
%
%   Input :
%           N : Number of samples
%           M : Number of measurements
%           K : Size of the support of x0
%           U and Ut : Synthesis/Analysis projection operators (e.g. Fourier, Hadamard, Noiselets ... etc)
%           Dyna : Dynamic range in dB
%           NoiseLevel : noise level in dB (i.e. SNR.  Higher means less
%               noise)
%           Mode : 
%               Mode = 0 will setup for k-sparse signals
%               Mode = 1 will setup for the first set of wavelet coefficents in approxsparse.mat
%               Mode = 2 will setup for the second set of wavelet coefficients
%           Repeat : a number (e.g. 1 - 10 ) that is unique for every test
%               this will influence the seed for the random number generator
%
%   Output :
%           x0 : N*1 double array - original signal
%           b : K*1 double array - data
%           A,At : Subsampling operator and its adjoint
%           supind : indices of the support of x0
%           Omega : indices of the observed projections
%           b_exact : the data, before noise has been added
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: February 2009
% Modified: Stephen Becker, Mon 3/30/09
% Modified: Stephen Becker, Nov 2009, adding support for complex data
%
% NESTA Version 1.1

function [x0,b,A,At,supind,Omega,b_exact]=msp_signal(N,M,K,U,Ut,Dyna,NoiseLevel,Mde,repeat)
global COMPLEX
if isempty(COMPLEX)
    COMPLEX = false;
end

if nargin < 8 || isempty(Mde)
    Mde = 0;  % Note: not Jerome's conventions
end

% Calculate the seed in a predetermined fashion
if Mde == 0
    seed = 500 + Dyna + repeat;
elseif Mde == 1
    seed = 650 + repeat;
else
    seed = 700 + repeat;
end

randn('state',seed); rand('state',seed);

if Mde
    supind=[];
    load approxsparse
    if Mde == 1, coeff = nat_wtc;
    elseif Mde == 2, coeff = nat_wtc2;
    else error('msp_signal: unknown Mde flag (0,1 or 2 are acceptable)'); end
    if N ~= length(coeff)
        disp('msp_signal: warning, N is not size of coefficients');
    end
    per = randperm( length(coeff) );
    x0 = coeff(per);
    x0 = x0(1:N);
else
    % the k-sparse case, with all synthetic data
    x0 = zeros(N,1);
    supind = randperm(N);
    supind = supind(1:K);
    valx = Dyna/20*rand(K,1);
    valx = valx - min(valx);valx = valx/max(valx)*Dyna/20;
    if ~COMPLEX
        x0(supind) = 10.^valx.*sign(randn(K,1));
    else
        x0(supind) = 10.^valx.*sign( randn(K,1) + 1i*randn(K,1) );
    end
end
Omega = randperm(N);
Omega = Omega(1:M);
A = @(x) x(Omega,:);
% A = @(z) A_Subset(z,Omega);
S.type = '()'; 
S.subs{1} = Omega; 
S.subs{2} = ':';
At = @(x) subsasgn( zeros(N,size(x,2)),S,x);
% At = @(z) At_Subset(z,N,Omega);
b = A(Ut(x0));
if nargout >= 7, b_exact = b; end
    
if ~COMPLEX
    b = b + NoiseLevel*randn(M,1);
else
    b = b + NoiseLevel*(randn(M,1)+1i*randn(M,1))/sqrt(2);
end
