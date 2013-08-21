%GMMB_COVFIXER   - force matrix to be a valid covariance matrix
%
% covmatrix = GMMB_COVFIXER(matrix)
% Matrix is forced (complex conjugate) symmetric,
% positive definite and its diagonal real valued.
%
% [covmatrix, loops] = GMMB_COVFIXER(...)
%    loops - number of rounds the positive definite fixer had to run.
%
% [covmatrix, loops, symerr] = GMMB_COVFIXER(...)
%    symerr - symmetry error matrix

%
% $Name:  $
% $Id: gmmb_covfixer.m,v 1.2 2004/11/02 09:00:18 paalanen Exp $
% Copyright 2003, Pekka Paalanen <pekka.paalanen@lut.fi>

% except isspd() function which is from The MathWorks Matlab mvnpdf.m.

function [nsigma, varargout] = gmmb_covfixer(sigma);

D = size(sigma, 1);
fixrate = 0.01;
covfixmat = ones(D) + fixrate*eye(D);
loops = 0;
min_limit = eps*10;

if ~all( isfinite(sigma(:))  )
	error('covariance matrix is not finite');
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
%ISPDS Test if a matrix is positive definite symmetric
% T = ISPDS(SIGMA) returns a logical indicating whether the matrix SIGMA is
% square, symmetric, and positive definite, i.e., it is a valid full rank
% covariance matrix.
%
% [T,R] = ISPDS(SIGMA) returns the cholesky factor of SIGMA in R.  If SIGMA
% is not square symmetric, ISPDS returns [] in R.

%   Copyright 1993-2002 The MathWorks, Inc.
%   Revision: 1.2   Date: 2002/03/28 16:51:27

% Test for square, symmetric
% NOTE: imagfixer already enforces squareness and symmetricity,
% and fixing affects only the diagonal, so this is not necessary
%[n,m] = size(Sigma);
%if (n == m) & all(all(abs(Sigma - Sigma') < 10*eps*max(abs(diag(Sigma)))));

    % Test for positive definiteness
    [R,p] = chol(Sigma);
    if p == 0
        t = 1;
    else
        t = 0;
    end

%else
%    R = [];
%    t = 0;
%end

% ------------------

function nsigma = imagfixer(sigma);

% force symmetric
nsigma = sigma - (sigma - sigma')/2;
% purge imag
purge = imag(diag(nsigma));
nsigma = nsigma - diag(purge)*1i;

if max(purge) > 1e-4
	warning_wrap('gmmbayes:covfixer:imagfixer', 'Quite big imaginary components removed from the diagonal');
end
