%GMMB_PDF    - (Complex range) multivariate Gaussian mixture model pdf
%
% p = gmmb_pdf(data, bayesS)
%
% data = N x D matrix
% bayesS = 1 x K struct array, the bayesS struct,
% fields used:
%   mu = D x C matrix
%   sigma = D x D x C matrix array
%   weight = C x 1 vector
%
% p = N x K matrix
%
% D dimensions, N points, C components, K classes
%
% The result is a matrix of PDF values computed for
% each data point and class.
% see: GMMB_CMVNPDF

% Author: Pekka Paalanen <pekka.paalanen@lut.fi>

%
% $Name:  $
% $Id: gmmb_pdf.m,v 1.2 2004/11/02 09:00:18 paalanen Exp $

function [pdf] = gmmb_pdf(data_, bayesS_);
N = size(data_,1);
K = size(bayesS_, 2);

pdf = zeros(N,K);

for k = 1:K
	weight = bayesS_(k).weight;
	mu = bayesS_(k).mu;
	sigma = bayesS_(k).sigma;

	p = zeros(N,1);
	for c = 1:size(weight,1);
		p = p + weight(c) * ...
		   gmmb_cmvnpdf(data_, mu(:,c).', sigma(:,:,c));
	end
	pdf(:,k) = p;
end

