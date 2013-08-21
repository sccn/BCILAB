% SB2_FULLSTATISTICS  Compute all relevant statistics in full for SPARSEBAYES
%
% [SIGMA,MU,S_IN,Q_IN,S_OUT,Q_OUT,FACTOR,LOGML,GAMMA,...
%	BETABASIS_PHI,BETA] = ...
%  SB2_FULLSTATISTICS(LIKELIHOOD,BASIS,PHI,TARGETS,USED,ALPHA,BETA,...
%       MU,BASIS_PHI,BASIS_TARGETS,OPTIONS)
%
% OUTPUT ARGUMENTS:
%
%	SIGMA			Posterior covariance matrix for relevant bases
%	MU				Posterior mean
%	S_IN			S-factors for in-model (relevant) basis vectors
%	Q_IN			Q-factors for in-model (relevant) basis vectors
%	S_OUT			S-factors for all basis vectors
%	Q_OUT			Q-factors for all basis vectors
%	FACTOR			Q^2-S relevance factors
%	LOGML			Log-marginal-likelihood
%	GAMMA			"Well-determinedness" factors
%	BETABASIS_PHI	Cached value of BASIS'*B*PHI matrix
%	BETA			Inverse noise variance (vector of beta
%					approximations in non-Gaussian case)
% 
% INPUT ARGUMENTS:
% 
%	LIKELIHOOD		LIKELIHOOD structure
%	BASIS			Full NxM basis matrix
%	PHI				Current relevant basis matrix
%	TARGETS			N-vector with target output values
%	USED			Relevant basis vector indices
%	ALPHA			Hyperparameters
%	BETA			Inverse noise variance (ignored in non-Gauss case)
%	MU				Current posterior mean (for non-Gauss)
%	BASIS_PHI		Cached value of BASIS'*PHI (for Gauss only)
%	BASIS_TARGETS	Cached value of BASIS'*TARGETS (for Gauss only)
%	OPTIONS			Standard OPTIONS structure (see SB2_USEROPTIONS)
% 
% NOTES:
%
% This function computes the posterior, and other, statistics for the
% SPARSEBAYES algorithm in "long-hand" fashion. i.e. it does not use the
% more effecient "iterative" updates.
% 
% It is required on every iteration (unless approximating) in the
% non-Gaussian case, and on iterations in the Gaussian case where the noise
% estimate (BETA) is updated.
% 
% This function is intended for internal use by SPARSEBAYES only.
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
function [SIGMA,Mu,S_in,Q_in,S_out,Q_out,Factor,logML,Gamma,...
	  betaBASIS_PHI,beta] = ...
    SB2_FullStatistics(LIKELIHOOD,BASIS,PHI,Targets,Used,Alpha,beta,...
		       Mu,BASIS_PHI,BASIS_Targets,OPTIONS)

% Mu is only required for non-Gauss
% BASIS_PHI and BASIS_Targets are only required for Gauss

% beta (a vector for non-Gauss) is returned only for non-Gauss

MAX_POSTMODE_ITS	= 25; % More than enough

[N M_full]	= size(BASIS);
M			= size(PHI,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE POSTERIOR

% 
% Compute full relevant statistics
% 
if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
  %
  % Posterior is analytic: use linear algebra here
  % 
  % Compute posterior covariance SIGMA (inverse Hessian) 
  % via Cholesky decomposition
  % 
  U		= chol(PHI'*PHI*beta + diag(Alpha));
  Ui	= inv(U);
  SIGMA	= Ui * Ui';
  % Posterior mean Mu
  Mu	= (SIGMA * (PHI'*Targets)) * beta;
  % Data error and likelihood
  y				= PHI * Mu;
  e				= (Targets - y);
  ED			= e'*e;
  %
  dataLikely	= (N*log(beta) - beta*ED)/2;
  %
else
  %
  % Posterior must be approximated: find the posterior mode as the basis
  % for the Laplace approximation
  % 
  [Mu U beta dataLikely] = ...
      SB2_PosteriorMode(LIKELIHOOD,PHI,Targets,Alpha,Mu, ...
						MAX_POSTMODE_ITS,OPTIONS);
  % Compute covariance approximation
  % 
  Ui	= inv(U);
  SIGMA	= Ui * Ui';
  % Compute posterior-mean-based outputs and errors
  % 
  switch LIKELIHOOD.InUse
   case LIKELIHOOD.Bernoulli,
    y	= SB2_Sigmoid(PHI * Mu);
   case LIKELIHOOD.Poisson
    y	= exp(PHI*Mu);
  end
  e		= (Targets-y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE LOG MARGINAL LIKELIHOOD

logdetHOver2	= sum(log(diag(U)));
logML			= dataLikely - (Mu.^2)'*Alpha/2 + ...
    sum(log(Alpha))/2 - logdetHOver2;
% Well-determinedness factors
DiagC	= sum(Ui.^2,2);
Gamma	= 1 - Alpha.*DiagC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE THE Q & S VALUES

%
% Q: "quality" factor - related to how well the basis function contributes
% to reducing the error
% 
% S: "sparsity factor" - related to how orthogonal a given basis function
% is to the currently used set of basis functions
% 
if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
  %
  % Gaussian case simple: beta a scalar
  % 
  betaBASIS_PHI	= beta*BASIS_PHI;			
  %
  % The S_in calculation exploits the fact that BASIS is normalised,
  %  i.e. sum(BASIS.^2,1)==ones(1,M_full)
  %
  S_in		= beta - sum((betaBASIS_PHI*Ui).^2,2);
  Q_in		= beta*(BASIS_Targets - BASIS_PHI*Mu);
else
  %
  % Non-Gaussian case: beta an N-vector
  % 
  betaBASIS_PHI	= BASIS'*(PHI.* (beta*ones(1,M)));
  S_in			= (beta'*(BASIS.^2))' - sum((betaBASIS_PHI*Ui).^2,2);
  Q_in			= BASIS'*e;
end
%
S_out		= S_in;
Q_out		= Q_in;
%
% S,Q with that basis excluded: equations (23)
% 
S_out(Used)	= (Alpha .* S_in(Used)) ./ (Alpha - S_in(Used));
Q_out(Used)	= (Alpha .* Q_in(Used)) ./ (Alpha - S_in(Used));
%
% Pre-compute the "relevance factor" for ongoing convenience
% 
Factor		= (Q_out.*Q_out - S_out);

