% data = GMMB_MKCPLX(mu, sigma, N)   Generate complex Gaussian data
%
% Generates N points of complex valued data according to
% complex Gaussian distribution with parameters (mu, sigma).
%
%  data  = N x D list of vectors (complex)
%  mu    = 1 x D vector (complex)
%  sigma = D x D matrix (complex, positive semi-definite)
%
% Author:
%   Pekka Paalanen <pekka.paalanen@lut.fi>
%
% References:
%   Y. L. Tong, The Multivariate Normal Distribution,
%    Springer Series in Statistics, Springer-Verlag New York, 1990
%    p. 185
%
% $Id: gmmb_mkcplx.m,v 1.1 2004/11/02 08:32:22 paalanen Exp $

function data = mkcplx(mu, sigma, N);

D = size(mu,2);

[R, p] = chol(sigma);
if p>0
	error('Sigma must be positive semi-definite.');
end

magn = randn(N, D);
phase = rand(N, D).*pi;

data = complex(magn.*cos(phase), magn.*sin(phase)) * R + repmat(mu, N, 1);
