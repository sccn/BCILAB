% SB2_INITIALISATION  Initialise everything for the SPARSEBAYES algorithm
%
% [LIKELIHOOD, BASIS, BASISSCALES, ALPHA, BETA, MU, PHI, USED] = ...
%    SB2_INITIALISATION(TYPE, BASIS, TARGETS, SETTINGS, OPTIONS)
%
% OUTPUT ARGUMENTS:
% 
%	LIKELIHOOD	Likelihood structure (from SB2_LIKELIHOODS)
%	BASIS		Pre-processed full basis matrix
%	BASISSCALES	Scaling factors from full basis pre-processing
%	ALPHA		Initial hyperparameter alpha values
%	BETA		Initial noise level (Gaussian)
%	MU			Initial weight values
%	PHI			Initial "relevant" basis
%	USED		Indices of initial basis vectors (columns of BASIS)
% 
% INPUT ARGUMENTS:
% 
%	TYPE		Likelihood: one of 'Gaussian', 'Bernoulli', 'Poisson'
%
%	BASIS		NxM matrix of all possible basis vectors 
%				(one column per basis function)
%
%	TARGETS		N-vector with target output values
% 
%	SETTINGS	Initialisation structure for main parameter values via
%				SB2_PARAMETERSETTINGS 
% 
%	OPTIONS		User options structure from SB2_USEROPTIONS
%	
% NOTES: 
% 
% This function initialises all necessary model parameters appropriately
% before entering the main loop of the SPARSEBAYES inference algorithm.
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
function [LIKELIHOOD, BASIS, BasisScales, Alpha, beta, Mu, PHI, Used] = ...
    SB2_Initialisation(likelihood_, BASIS, Targets, SETTINGS, OPTIONS)

% A "reasonable" initial value for the noise in the Gaussian case
GAUSSIAN_SNR_INIT	= 0.1;

% "Reasonable" initial alpha bounds
INIT_ALPHA_MAX	= 1e3;
INIT_ALPHA_MIN	= 1e-3;

% 
% BASIS PREPROCESSING: 
% 
% Scale basis vectors to unit norm. This eases some calculations and 
% will improve numerical robustness later.
% 
[BASIS, BasisScales] = SB2_PreProcessBasis(BASIS);

%% Noise model considerations
 
% Validate the likelihood model specification
% 
LIKELIHOOD	= SB2_Likelihoods(likelihood_);

% Default beta value for non-Gaussian likelihood
beta		= [];
%
% In the Gaussian case, initialise "sensibly" to a proportion of the
% signal level (e.g. 10 percent)
% 
if LIKELIHOOD.InUse == LIKELIHOOD.Gaussian
  % Noise initialisation
  if ~isempty(SETTINGS.beta)
    % The user set "beta"
    beta	= SETTINGS.beta;
  elseif ~isempty(SETTINGS.noiseStdDev) && SETTINGS.noiseStdDev>0
    % The user set the "noise std dev"
    beta	= 1/SETTINGS.noiseStdDev^2;
  else
    % catch the pathological case where all outputs are zero
	% (although we're probably doomed anyway if that's true)
    stdt	= max([1e-6 std(Targets)]);
    % Initialise automatically approximately according to "SNR"
    beta	= 1/(stdt*GAUSSIAN_SNR_INIT)^2;
  end
end


% Initialise basis (PHI), mu and alpha
%
% Either as specified by the SETTINGS structure, or heuristically
% 

% First, compute 'linearised' output for use in heuristic initialisation
%
TargetsPseudoLinear	= Targets; % standard linear case
switch LIKELIHOOD.InUse
 case LIKELIHOOD.Bernoulli,
  TargetsPseudoLinear	= (2*Targets-1);
 case LIKELIHOOD.Poisson,
  TargetsPseudoLinear	= log(Targets+1e-3);
end


% 1) the starting basis, PHI

% Take account of "free basis": it needs to be included from the outset
% 
Extra	= setdiff(OPTIONS.freeBasis, SETTINGS.Relevant);
%
% [TODO:] More elegant to put "Extra" first in Used list?
% 
Used	= [SETTINGS.Relevant Extra];
%
% At this point Used will contain both pre-specified relevant basis
% functions (SB2_PARAMETERSETTINGS) and any "free" ones (SB2_USEROPTIONS).
% If neither are active, Used will be empty.
% 
if isempty(Used)
  % Set initial basis to be the largest projection with the targets
  proj			= (BASIS'*TargetsPseudoLinear);
  [foo Used]	= max(abs(proj));
  SB2_Diagnostic(OPTIONS, 2, ...
				 ['Initialising with maximally aligned ' ...
				  'basis vector (%d)\n'], Used);
else  
  SB2_Diagnostic(OPTIONS, 2, ...
				 'Initialising with supplied basis of size M = %d\n', ...
				 length(Used));
end
PHI	= BASIS(:,Used);
M	= length(Used);

% 2) the most probable weights

if isempty(SETTINGS.Mu)
  switch LIKELIHOOD.InUse
    %
   case LIKELIHOOD.Gaussian,
    % Mu will be analytically calculated later in the Gaussian case
    %
    Mu		= [];
   case LIKELIHOOD.Bernoulli,
    % Heuristic initialisation based on log-odds
    LogOut	= (TargetsPseudoLinear*0.9+1)/2;
    Mu		= PHI \ (log(LogOut./(1-LogOut)));
    %
   case LIKELIHOOD.Poisson
    % Heuristic initialisation based on log
    Mu		= PHI \ TargetsPseudoLinear;
    %
  end
else
  if length(SETTINGS.Mu)~=length(SETTINGS.Relevant)
    error('Basis length (%d) should equal weight vector length (%d)', ...
		  length(SETTINGS.Relevant), length(SETTINGS.Mu))
  end
  SB2_Diagnostic(OPTIONS, 2,'Initialising with supplied weights\n');
  % We may have extra "freebasis" functions without supplied weights
  % - set those to zero for now
  Mu	= [Mu ; zeros(length(Extra),1)];
end

% 3) the hyperparameters

if isempty(SETTINGS.Alpha)
  switch LIKELIHOOD.InUse
    %
   case LIKELIHOOD.Gaussian,
    % Exact for single basis function case (diag irrelevant), 
    % heuristic in the multiple case
    p		= diag(PHI'*PHI)*beta;
    q		= (PHI'*Targets)*beta;
    Alpha	= p.^2./(q.^2-p); %
    % Its possible that one of the basis functions could already be
    % irrelevant (alpha<0), so trap that
    if all(Alpha<0)
      SB2_Diagnostic(OPTIONS, 1,...
		     ['WARNING: no relevant basis function at' ...
		      ' initialisation\n']);
    end
	% The main algorithm will handle these automatically shortly
	% (i.e. prune them)
    Alpha(Alpha<0) = INIT_ALPHA_MAX; 
    %
   case {LIKELIHOOD.Bernoulli, LIKELIHOOD.Poisson}
    % Heuristic initialisation, trapping Mu==0
    Alpha						= 1./(Mu+(Mu==0)).^2;
    % Limit max/min alpha
    Alpha(Alpha<INIT_ALPHA_MIN)	= INIT_ALPHA_MIN;
    Alpha(Alpha>INIT_ALPHA_MAX)	= INIT_ALPHA_MAX;
    %
  end
  if M==1
    SB2_Diagnostic(OPTIONS,2,'Initial alpha = %g\n', Alpha);
  end
else
  if length(SETTINGS.Alpha)~=length(SETTINGS.Relevant)
    error('Basis length (%d) should equal alpha vector length (%d)', ...
		  length(SETTINGS.Relevant), length(SETTINGS.Alpha))
  end
  % We may have extra "freebasis" functions without supplied alpha
  % - set those to zero for now
  Alpha	= [SETTINGS.Alpha ; zeros(length(Extra),1)];
end

% Finally - set Alphas to zero for free basis
ALPHA_ZERO								= eps;
Alpha(ismember(Used,OPTIONS.freeBasis))	= ALPHA_ZERO;