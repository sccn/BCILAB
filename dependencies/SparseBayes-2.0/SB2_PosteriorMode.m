% SB2_POSTERIORMODE  Posterior mode-finder for the SPARSEBAYES algorithm
%
% [MU, U, BETA, LIKELIHOODMODE, BADHESS] = ...
%    SB2_POSTERIORMODE(LIKELIHOOD,BASIS,TARGETS,ALPHA,MU,ITSMAX,OPTIONS)
%
% OUTPUT ARGUMENTS:
% 
%	MU				Parameter values at the mode
%	U				Cholesky factor of the covariance at the mode
%	BETA			Vector of pseudo-noise variances at the mode
%	LIKELIHOODMODE	Data likelihood at the mode
%	BADHESS			Returns true if Hessian is "bad" (becomes
%					non-positive-definite during maximisation)
% 
% INPUT ARGUMENTS:
% 
%	LIKELIHOOD	LIKELIHOOD structure
%	BASIS		Current relevant basis matrix
%	TARGETS		N-vector with target output values
%	ALPHA		Current hyperparameters
%	MU			Current weights
%	ITSMAX		Maximum number of iterations to run
%	OPTIONS		Standard OPTIONS structure (only affects diagnostics)
% 
% NOTES:
% 
% SB2_POSTERIORMODE finds the posterior mode (with respect to the
% weights) of the likelihood function in the non-Gaussian case to
% facilitate subsequent Laplace approximation.
%
% This function is intended for internal use by SPARSEBAYES only (within
% SB2_FullStatistics).
%

%
% Copyright 2009, Vector Anomaly Ltd
%
% This file is part of the SPARSEBAYES library for Matlab (V2.0).
%
% SPARSEBAYES is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
%
% SPARSEBAYES is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with SPARSEBAYES in the accompanying file "licence.txt"; if not, write to
% the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
% MA 02110-1301 USA
%
% Contact the author: m a i l [at] m i k e t i p p i n g . c o m
%
function [Mu, U, Beta, likelihoodMode, badHess] = ...
    SB2_PosteriorMode(LIKELIHOOD,BASIS,Targets,Alpha,Mu,itsMax,OPTIONS)

% TOLERANCES
% 
% Termination criterion for each gradient dimension
% 
GRADIENT_MIN	= 1e-6;
%
% Minimum fraction of the full Newton step considered
%  
STEP_MIN	= 1/(2^8);
%

[N, M]	= size(BASIS);
A		= diag(Alpha);

% NB: for historical reasons, we work in term of error here (negative
% log-liklihood) and minimise

%
% Get current model output and data error
% 
BASIS_Mu		= BASIS*Mu; % Linear output
[dataError, y]	= SB2_DataError(LIKELIHOOD, BASIS_Mu, Targets);
%
% Add on the weight penalty
%
regulariser		= (Alpha'*(Mu.^2))/2;
newTotalError	= dataError + regulariser;
%
badHess		= false;
errorLog	= zeros(itsMax,1);
%
for iteration=1:itsMax
  %
  % Log the error value each iteration
  %
  errorLog(iteration)	= newTotalError;
  SB2_Diagnostic(OPTIONS,4,'PM cycle: %2d\t error: %.6f\n', ...
				 iteration, errorLog(iteration));

  % Construct the gradient
  %
  e	= (Targets-y);
  g	= BASIS'*e - Alpha.*Mu;
  %
  % Compute the likelihood-dependent analogue of the noise precision.
  % NB: Beta now a vector.
  % 
  switch LIKELIHOOD.InUse
   case LIKELIHOOD.Bernoulli
    Beta	= y.*(1-y);
   case LIKELIHOOD.Poisson
    Beta	= y;
  end
  % Compute the Hessian
  BASIS_B	= BASIS .* (Beta * ones(1,M));
  H			= (BASIS_B'*BASIS + A);
  % Invert Hessian via Cholesky, watching out for ill-conditioning
  [U, pdErr]	= chol(H);
  % Make sure its positive definite
  if pdErr
	% If you see this, it's *probably* the result of a bad choice of
    % basis. e.g. a kernel with too large a "width"
	% 
    SB2_Diagnostic(OPTIONS, 1,...
				   '** Warning ** Ill-conditioned Hessian (%g)\n', ...
				   rcond(H));
    badHess			= true;
    U				= [];
    Beta			= [];
    likelihoodMode	= [];
    return
  end
  %
  % Before progressing, check for termination based on the gradient norm
  % 
  if all(abs(g)<GRADIENT_MIN)
    errorLog	= errorLog(1:iteration);
    SB2_Diagnostic(OPTIONS,4,['PM convergence (<%g) after ' ...
					'%d iterations, |g| = %g\n'], ...
				   GRADIENT_MIN,iteration,max(abs(g)));
    break
  end

  % If all OK, compute full Newton step: H^{-1} * g
  % 
  DeltaMu	= U \ (U' \ g);
  step		= 1;
  %
  while step>STEP_MIN
	% Follow gradient to get new value of parameters
	% 
    Mu_new		= Mu + step*DeltaMu;
    BASIS_Mu	= BASIS*Mu_new;
	%
    % Compute outputs and error at new point
	% 
    [dataError,y]	= SB2_DataError(LIKELIHOOD,BASIS_Mu,Targets);
    regulariser		= (Alpha'*(Mu_new.^2))/2;
    newTotalError	= dataError + regulariser;
	%
    % Test that we haven't made things worse
	% 
    if newTotalError>=errorLog(iteration)
      % If so, back off!
	  % 
      step	= step/2;
      SB2_Diagnostic(OPTIONS, 4,['PM error increase! Backing off to l=' ...
					'%.3f\n'], step);
    else
      Mu	= Mu_new;
      step	= 0;			% this will force exit from the "while" loop
    end
  end
  %
  % If we get here with non-zero "step", it means that the smallest
  % offset from the current point along the "downhill" direction did not
  % lead to a decrease in error. In other words, we must be
  % infinitesimally close to a minimum (which is OK).
  % 
  if step
    SB2_Diagnostic(OPTIONS, 4, ...
				   'PM stopping due to back-off limit (|g| = %g)\n', ...
				   max(abs(g)));
    break;
  end
end
%
% Simple computation of return value of log likelihood at mode
% 
likelihoodMode	= -dataError;

%%%%%%%%%%%%%%%%%%%%%%
% 
% Support function
%
%%%%%%%%%%%%%%%%%%%%%%

function [e,y] = SB2_DataError(LIKELIHOOD,BASIS_Mu,Targets)

switch LIKELIHOOD.InUse
 case LIKELIHOOD.Bernoulli
  y	= SB2_Sigmoid(BASIS_Mu);
  % Handle probability zero cases
  y0	= (y==0);
  y1	= (y==1);
  if any(y0(Targets>0)) || any(y1(Targets<1))
    % Error infinite when model gives zero probability in
    % contradiction to data
    e	= inf;
  else
    % Any y=0 or y=1 cases must now be accompanied by appropriate
    % output=0 or output=1 values, so can be excluded.
    e	= -(Targets(~y0)'*log(y(~y0)) + (1-Targets(~y1))'*log(1-y(~y1)));
  end
  %
 case LIKELIHOOD.Poisson
  y	= exp(BASIS_Mu);
  e	= -sum(Targets.*BASIS_Mu - y);
  %
end
