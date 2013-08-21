%GMMB_CMVNPDF - Compute the value of Gaussian PDF (real or complex range)
%
%     Y = GMMB_CMVNPDF(X, MU, SIGMA)
%     Computes the D-dimensional (complex) Gaussian PDF with parameters
%     MU and SIGMA in points X(i,:) -> Y(i), i=1..N.
%     X:     N x D matrix of row vectors
%     MU:    1 x D vector
%     SIGMA: D x D matrix
%     SIGMA is assumed complex conjugate symmetric and semi-positive definite.
%        (You can use CHOL to test it, or gmmb_covfixer.)
%
%    This is a dual function with different formulas whether MU is
%    complex variable or not.
%
%    NOTE: This code is not based on the Matlab(r) mvnpdf.m and so
%    does not include the same sanity checks.
%
% Author(s):
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%
% References:
%   [1] Albertazzi, G., Cioni, S., Corazza, G.E., Vanelli-Coralli, A.,
%    Turbo Code Performance over Fading Channels in S-UMTS,
%    Research Center on Advanced Electronic Systems for Information
%    and Communication Technologies, University of Bologna
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf is
%   Copyright (C) 2003, 2004 by Pekka Paalanen and Joni-Kristian
%   Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2004/11/02 09:00:18 $
%

function y = gmmb_cmvnpdf(X, Mu, Sigma)

% Get size of data.
[n,d] = size(X);

invSigma = inv(Sigma);

invSigmaMu = invSigma*Mu';

sumvec = sum((X*invSigma).*conj(X),2);

sqrdist = sumvec ...
	- X*invSigmaMu ...
	- ((Mu*invSigma)*X').' ...
	+ Mu*invSigmaMu;

invDetSigma = 1/real(det(Sigma));

if isreal(Mu)
	y = sqrt( (2*pi)^(-d) * invDetSigma ) .* exp(-0.5*real(sqrdist));
else
	y = ( pi^(-d) * invDetSigma ) .* exp(-real(sqrdist));
end
