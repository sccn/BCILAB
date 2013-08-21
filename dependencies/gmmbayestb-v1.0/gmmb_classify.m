%GMMB_CLASSIFY Classify data using Bayesian or Mahalanobis distance classifier.
%
%     T = GMMB_CLASSIFY(S, data, ...) Classifies D dimensional data (N points)
%     using Gaussian Mixture Model
%     Bayesian classifier in struct S into K classes.
%     S is a bayesS struct, see readme.txt.
%
%     See also GMMB_CREATE.
%
%     ***
%     This is a legacy interface that is no longer developed.
%     This does not use gmmb_generatehist but relies on the histS
%     created during classifier training.
%     ***
%
% Optional parameters:
%    values        Return numbers instead of class labels.
%                  For Mahalanobis boolean; return the Mahalanobis distances.
%                  For Bayesian: 0 class labels
%                                1 posterior probabilities
%                                2 posterior likelihood (no scaling)
%                  Default setting 0. T: N x K matrix, if enabled.
%
%    mahalanobis   Use Mahalanobis distances as classification rule
%                  instead of Bayesian. Boolean, default 0.
%
%    quantile      Use a 'trash class' based on training data density quantile.
%                  The function gmmb_generatehist will be called.
%                  Not available with Mahalanobis distance method.
%                  Range [0, 1], default 0 (not in use).
%                  If a sample goes to the trash class:
%                    class label = 0
%                    all posterior probabilities = 0
%                    all posterior likelihoods = 0
%                  If a sample does not belong to a class in case values~=0,
%                  the respective posterior will be zero.
%
% References:
%   [1] Duda, R.O., Hart, P.E, Stork, D.G, Pattern Classification,
%   2nd ed., John Wiley & Sons, Inc., 2001.
%
% Author(s):
%    Joni Kamarainen <Joni.Kamarainen@lut.fi>
%    Pekka Paalanen <pekka.paalanen@lut.fi>
%
% Copyright:
%
%   Bayesian Classifier with Gaussian Mixture Model Pdf
%   functionality is Copyright (C) 2003, 2004 by Pekka Paalanen and
%   Joni-Kristian Kamarainen.
%
%   $Name:  $ $Revision: 1.2 $  $Date: 2005/04/14 10:33:34 $
%

function [t] = gmmb_classify(bayesS, data_, varargin);

conf = struct(...
	'values', 0, ...
	'mahalanobis', 0, ...
	'quantile', 0);
conf = getargs(conf, varargin);


N = size(data_,1);
K = size(bayesS,2);

% data_ is N x D

if conf.quantile ~= 0
	histS = gmmb_generatehist(bayesS, 1000);
	conf_thresh = gmmb_frac2lhood(histS, ...
	                  conf.quantile*ones(1,K));
	clear histS;
end


if conf.mahalanobis
	%disp('Using Mahalanobis distance.');
	% Mahalanobis distance classifier
	sqrmdist = zeros(N,K);
	for k = 1:K
		C = size(bayesS(k).mu, 2);
		sqrdist = zeros(N,C);
		for c = 1:C
			invs = inv(bayesS(k).sigma(:,:,c));
			mu = bayesS(k).mu(:,c).';
			sqrdist(:,c) = sum((data_*invs).*conj(data_),2) ...
				- data_*invs*mu' ...
				- (mu*invs*data_').' ...
				+ mu*invs*mu';
		end
		sqrmdist(:,k) = min(real(sqrdist), [], 2);
	end
	if conf.values ~= 0
		t = sqrmdist;
	else
		[a, b] = min(sqrmdist, [], 2);
		t = b';
	end
else
	% GMM Bayesian classifier
	
	% classify all points simultaneously
	pxomega = gmmb_pdf(data_, bayesS);
	tulo = gmmb_weightprior(pxomega, bayesS);
	% tulo is the product of GMM PDF values and class apriories
	
	% Zero out class likelihoods of samples that do not belong
	% to the class, based on the density quantile.
	if conf.quantile ~= 0
		mask = (pxomega < repmat(conf_thresh, N, 1));
		tulo(mask) = 0;
	end
	
	% Compute posteriors if requested.
	if conf.values == 1
		t = gmmb_normalize(tulo);
	else
		t = tulo;
	end

	% Find the classification outcome if requested
	% Mark the samples that got into none of the known classes
	if (conf.values == 0)
		t = gmmb_decide(t);
	end
end

