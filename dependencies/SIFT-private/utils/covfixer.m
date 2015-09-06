%COVFIXER   - force matrix to be a valid covariance matrix
%
% covmatrix = COVFIXER(matrix)
% Matrix is forced (complex conjugate) symmetric,
% positive definite and its diagonal real valued.
%
% [covmatrix, loops] = COVFIXER(...)
%    loops - number of rounds the positive definite fixer had to run.
%
% [covmatrix, loops, symerr] = COVFIXER(...)
%    symerr - symmetry error matrix

%
% $Name:  $
% $Id: covfixer.m,v 1.3 2005/06/16 06:12:14 paalanen Exp $
% Copyright 2003, Pekka Paalanen <pekka.paalanen@lut.fi>

function [nsigma, varargout] = covfixer(sigma)

D = size(sigma, 1);
fixrate = 0.01;
covfixmat = ones(D) + fixrate*eye(D);
loops = 0;
min_limit = eps*10;

if ~all( isfinite(sigma(:))  )
	warning('gmmbayes:covfixer','covariance matrix is not finite');
end

% Running imagfixer is not counted as covariance fixing,
% the changes are assumed to be so small.
nsigma = imagfixer(sigma);

if nargout>2
	varargout(2) = {(sigma-nsigma)};
end

while isspd(nsigma) == 0
	% covariance matrix is not positive definite
	% fix it
	loops  = loops+1;
	d = diag(nsigma);
	if any(d <= min_limit)
		% negative or zero (<eps) on the diagonal
		m = max(abs(d)) * fixrate;
		neg = min(d);
		if neg < 0
			% there is a negative component on the diagonal
			% get rid of it.
			addit = (m-neg)*eye(D);
		else
			if m < min_limit
				m = min_limit;
			end
			addit = m*eye(D);
		end
		nsigma = nsigma + addit;
	else
		% increase diagonal values by 1 percent
		nsigma = nsigma .* covfixmat;
	end
end

if nargout>1
	varargout(1) = {loops};
end


% ------------------

function [t,R] = isspd(Sigma)
%ISSPD Test if a matrix is positive definite symmetric

% Test for positive definiteness
t = 1;

[R,p] = chol(Sigma);
if p ~= 0
	t = 0;
	return;
end

if any( svd(Sigma) < 10*eps )
	t = 0;
end

% ------------------

function nsigma = imagfixer(sigma)

% force symmetric
nsigma = sigma - (sigma - sigma')/2;
% purge imag
purge = imag(diag(nsigma));
nsigma = nsigma - diag(purge)*1i;

if max(purge) > 1e-4
	warning_wrap('gmmbayes:covfixer:imagfixer', 'Quite big imaginary components removed from the diagonal');
end