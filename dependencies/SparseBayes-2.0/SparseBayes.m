% SPARSEBAYES  Sparse Bayesian modelling: main estimation algorithm
%
% [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = ...
%    SPARSEBAYES(LIKELIHOOD, BASIS, TARGETS, OPTIONS, SETTINGS)
%
% OUTPUT ARGUMENTS:
% 
%	PARAMETER	Structure specifying inferred primary parameters:
% 
%	.Value		Vector of weight values
%	.Relevant	Vector of corresponding indices of relevant
%				columns of BASIS matrix (sorted ascending)
%
%	HYPERPARAMETER	Structure specifying inferred hyperparameters:
% 
%	.Alpha			Vector of weight precision values
%	.beta			Noise precision (Gaussian likelihood case)
% 
%	DIAGNOSTIC	Structure containing various diagnostics:
% 
%	.Gamma		Vector of "well-determined" factors [0,1] for
%				relevant weights 
%	.Likelihood	Vector of evolving log-marginal-likelihood
%	.iterations	Number of iterations run
%	.S_Factor	Vector of S ("Sparsity") values for relevant weights
%	.Q_Factor	Vector of Q ("Quality") values for relevant weights
% 
% INPUT ARGUMENTS:
% 
%	LIKELIHOOD	String comprising one of 'Gaussian', 'Bernoulli' or 'Poisson'
%
%	BASIS		NxM matrix of basis vectors (one column per basis function)
%
%	TARGETS		N-vector with target output values
% 
%	OPTIONS		User-specified settings via SB2_USEROPTIONS [Optional]
% 
%	SETTINGS	Initialisation of main parameter values via
%				SB2_PARAMETERSETTINGS [Optional]
%
% NOTES: 
%
% SPARSEBAYES is the implementation of the main algorithm for parameter
% inference in "sparse Bayesian" models.
% 
% Given inputs (BASIS), desired outputs (TARGETS) and an appropriate
% LIKELIHOOD function, SPARSEBAYES will optimise the log-marginal-likelihood
% of the corresponding sparse Bayesian model and should return (given
% reasonable choice of basis) a sparse vector of model parameters.
% 
% OPTIONS and SETTINGS arguments may be omitted, and if so, will assume
% sensible default values.
% 
% SEE ALSO:
%
%	SB2_USEROPTIONS: available user-definable options, such as the
%	number of iterations to run for, level of diagnostic output etc.
% 
%	SB2_PARAMETERSETTINGS: facility to change default initial values for
%	parameters and hyperparameters.
% 
%	SB2_CONTROLSETTINGS: hard-wired internal algorithm settings.
%
% The main algorithm is based upon that outlined in "Fast marginal
% likelihood maximisation for sparse Bayesian models", by Tipping & Faul, in
% Proceedings of AISTATS 2003. That paper may be downloaded from
% www.relevancevector.com or via the conference online proceedings site at
% http://research.microsoft.com/conferences/aistats2003/proceedings/.
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
function [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = ...
    SparseBayes(likelihood_, BASIS, Targets, OPTIONS, SETTINGS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SET UP INITIAL PARAMETERS, USER OPTIONS AND ALGORITHM CONTROL SETTINGS
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If no initial parameter setting structure passed, import the defaults
if ~exist('SETTINGS','var')
  SETTINGS	= SB2_ParameterSettings;
end

%% If no user options passed, import defaults
if ~exist('OPTIONS','var')
  OPTIONS	= SB2_UserOptions;
end

%% Any sanity checks on options and initialisation here

% Error if fixed noise specified but value not set
% 
if OPTIONS.fixedNoise && isempty(SETTINGS.beta) && ...
    isempty(SETTINGS.noiseStdDev)
  error('Option to fix noise variance is set but value is not supplied.')
end

%% Get the default algorithm control settings

CONTROLS	= SB2_ControlSettings;

% Start the clock now for diagnostic purposes (and to allow the algorithm to
% run for a fixed time)
% 
t0			= clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% INITIALISATION
%%
%% Pre-process basis, set up internal parameters according to initial
%% settings and likelihod specification
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kick off diagnostics (primarily, open log file if specified)
% 
OPTIONS = SB2_Diagnostic(OPTIONS, 'start');

%
% Initialise everything, based on SETTINGS and OPTIONS
% 
[LIKELIHOOD, BASIS, BasisScales, Alpha, beta, Mu, PHI, Used] = ...
    SB2_Initialisation(likelihood_, BASIS, Targets, SETTINGS, OPTIONS);
%
% Cache some values for later efficiency
% 
if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
  % It will be computationally advantageous to "cache" this quantity 
  % in the Gaussian case
  BASIS_PHI	= BASIS'*PHI;
else
  BASIS_PHI	= [];
end
BASIS_Targets	= BASIS'*Targets;

% FULL COMPUTATION
% 
% Initialise with a full explicit computation of the statistics
% 
% NOTE: The AISTATS paper uses "S/Q" (upper case) to denote the key
% "sparsity/quality" Factors for "included" basis functions, and "s/q"
% (lower case) for the factors calculated when the relevant basis
% functions are "excluded".
% 
% Here, for greater clarity:
% 
%	All S/Q are denoted by vectors S_in, Q_in
%	All s/q are denoted by vectors S_out, Q_out
% 
[SIGMA,Mu,S_in,Q_in,S_out,Q_out,Factor,logML,Gamma,BASIS_B_PHI,beta] = ...
      SB2_FullStatistics(LIKELIHOOD,BASIS,PHI,Targets,Used, ...
			 Alpha,beta,Mu,BASIS_PHI,BASIS_Targets,OPTIONS); 

%
% Avoid falling over in pathological case of zero iterations
% 
if OPTIONS.iterations==0
  PARAMETER				= [];
  HYPERPARAMETER		= [];
  DIAGNOSTIC.Likelihood = logML; 
  return
end

%
[N M_full]	= size(BASIS);
M			= size(PHI,2);
%
% Some diagnostics
% 
addCount	= 0;
deleteCount	= 0;
updateCount	= 0;
%
% Create storage to record evolution of log marginal likelihood
% 
maxLogSize		= OPTIONS.iterations + CONTROLS.BetaUpdateStart + ...
    ceil(OPTIONS.iterations/CONTROLS.BetaUpdateFrequency);
logMarginalLog	= zeros(maxLogSize,1);
count			= 0;
%
% If we're doing basis alignment testing, we'll need to maintain lists of
% those functions that are near identical, both in and out of the current
% model.
% 
if CONTROLS.BasisAlignmentTest
  Aligned_out		= [];
  Aligned_in		= [];
  alignDeferCount	= 0;
end

% ACTION CODES
%
% Assign an integer code to the basic action types
% 
ACTION_REESTIMATE	= 0;			
ACTION_ADD			= 1;
ACTION_DELETE		= -1;
%
% Some extra types
%
ACTION_TERMINATE		= 10;
%
ACTION_NOISE_ONLY		= 11;
%
ACTION_ALIGNMENT_SKIP	= 12;

%
% Before kicking off the main loop, call the specified "callback" function
% with "ACTION_ADD" to take account of the initialisation of the model
%
if OPTIONS.callback
  feval(OPTIONS.callbackFunc, 0, ACTION_ADD,logML/N, ...
		Used, Mu./BasisScales(Used)', SIGMA,...
		Alpha, beta, Gamma, PHI,...
		OPTIONS.callbackData{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% MAIN LOOP
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
i				= 0;	% Iteration number
LAST_ITERATION	= false;
%
while (~LAST_ITERATION)
  %
  i	= i+1;
  %
  % "UpdateIteration": set to true if this is an iteration where fast matrix
  % update rules can be used compute the appropriate quantities
  % 
  % This can be done if:
  % 
  % -	we are using a Gaussian likelihood
  % -	we are using other likelihoods and we have not specified a full
  %		posterior mode computation this cycle
  % 
  UpdateIteration	= LIKELIHOOD.InUse==LIKELIHOOD.Gaussian || ...
      rem(i,CONTROLS.PosteriorModeFrequency);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% DECISION PHASE
  %%
  %% Assess all potential actions
  %%
  
  %
  % Compute change in likelihood for all possible updates
  % 
  DeltaML		= zeros(M_full,1);
  %
  Action		= ACTION_REESTIMATE*ones(M_full,1); % Default
  %
  % 'Relevance Factor' (Q^S-S) values for basis functions in model
  % 
  UsedFactor	= Factor(Used);
  
  %
  % RE-ESTIMATION: must be a POSITIVE 'factor' and already IN the model
  % 
  iu		= UsedFactor>CONTROLS.ZeroFactor;
  index		= Used(iu);
  NewAlpha	= S_out(index).^2 ./ Factor(index);
  Delta		= (1./NewAlpha - 1./Alpha(iu)); % Temp vector
  %
  % Quick computation of change in log-likelihood given all re-estimations
  % 
  DeltaML(index)	= (Delta.*(Q_in(index).^2) ./ ...
			   (Delta.*S_in(index) + 1) - ...
			   log(1 + S_in(index).*Delta))/2;
  
  % 
  % DELETION: if NEGATIVE factor and IN model
  %
  % But don't delete:
  %		- any "free" basis functions (e.g. the "bias")
  %		- if there is only one basis function (M=1)
  % 
  % (In practice, this latter event ought only to happen with the Gaussian
  % likelihood when initial noise is too high. In that case, a later beta
  % update should 'cure' this.)
  % 
  iu			= ~iu; 	% iu = UsedFactor <= CONTROLS.ZeroFactor
  index			= Used(iu);
  anyToDelete	= ~isempty(setdiff(index,OPTIONS.freeBasis)) && M>1;
  %
  if anyToDelete
	%
	% Quick computation of change in log-likelihood given all deletions
	% 
    DeltaML(index)	= -(Q_out(index).^2 ./ (S_out(index) + Alpha(iu)) - ...
			    log(1 + S_out(index) ./ Alpha(iu)))/2;
    Action(index)	= ACTION_DELETE;
    % Note: if M==1, DeltaML will be left as zero, which is fine
  end
  
  % 
  % ADDITION: must be a POSITIVE factor and OUT of the model
  % 
  % Find ALL good factors ...
  GoodFactor		= Factor>CONTROLS.ZeroFactor;
  % ... then mask out those already in model
  GoodFactor(Used)	= 0;		
  % ... and then mask out any that are aligned with those in the model
  if CONTROLS.BasisAlignmentTest
    GoodFactor(Aligned_out)	= 0;
  end
  %
  index			= find(GoodFactor);
  anyToAdd		= ~isempty(index);
  if anyToAdd
	%
	% Quick computation of change in log-likelihood given all additions
	% 
    quot			= Q_in(index).^2 ./ S_in(index);
    DeltaML(index)	= (quot - 1 - log(quot))/2;
    Action(index)	= ACTION_ADD;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Post-process action results to take account of preferences
  
  % Ensure that nothing happens with "free basis" functions
  % 
  DeltaML(OPTIONS.freeBasis)	= 0;
  
  % If we prefer ADD or DELETE actions over RE-ESTIMATION
  % 
  if (anyToAdd && CONTROLS.PriorityAddition) || ...
      (anyToDelete && CONTROLS.PriorityDeletion)
    % We won't perform re-estimation this iteration, which we achieve by
    % zero-ing out the delta
    DeltaML(Action==ACTION_REESTIMATE)	= 0;
    % Furthermore, we should enforce ADD if preferred and DELETE is not
    % - and vice-versa
    if (anyToAdd && CONTROLS.PriorityAddition && ~CONTROLS.PriorityDeletion)
      DeltaML(Action==ACTION_DELETE)	= 0;
    end
    if (anyToDelete && CONTROLS.PriorityDeletion && ~CONTROLS.PriorityAddition)
      DeltaML(Action==ACTION_ADD)		= 0;
    end
  end
  
  % Finally...we choose the action that results 
  % in the greatest change in likelihood
  % 
  [deltaLogMarginal nu]	= max(DeltaML);
  selectedAction		= Action(nu);
  anyWorthwhileAction	= deltaLogMarginal>0;
  %
  % We need to note if basis nu is already in the model, and if so,
  % find its interior index, denoted by "j"
  %
  if selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE	
    j		= find(Used==nu);
  end
  %
  % Get the individual basis vector for update and compute its optimal alpha,
  % according to equation (20): alpha = S_out^2 / (Q_out^2 - S_out) 
  %
  Phi		= BASIS(:,nu);
  newAlpha	= S_out(nu)^2 / Factor(nu);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % TERMINATION CONDITIONS
  %
  % Propose to terminate if:
  % 
  % 1.	there is no worthwhile (likelihood-increasing) action, OR
  % 
  % 2a.	the best action is an ACTION_REESTIMATE but this would only lead to
  %		an infinitesimal alpha change, AND
  % 2b.	at the same time there are no potential awaiting deletions
  % 
  if ~anyWorthwhileAction || ...
	(selectedAction==ACTION_REESTIMATE && ...
	 abs(log(newAlpha) - log(Alpha(j)))<CONTROLS.MinDeltaLogAlpha && ...
	 ~anyToDelete)
	%
	selectedAction	= ACTION_TERMINATE;
	act_			= 'potential termination';
	%
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % ALIGNMENT CHECKS
  %
  % If we're checking "alignment", we may have further processing to do
  % on addition and deletion
  % 
  if CONTROLS.BasisAlignmentTest
    %
    % Addition - rule out addition (from now onwards) if the new basis
    % vector is aligned too closely to one or more already in the model
    % 
    if selectedAction==ACTION_ADD
      % Basic test for correlated basis vectors
	  % (note, Phi and columns of PHI are normalised)
	  % 
      p				= Phi'*PHI;
      findAligned	= find(p>CONTROLS.AlignmentMax);
      numAligned	= length(findAligned);
      if numAligned>0
		% The added basis function is effectively indistinguishable from
        % one present already
		selectedAction	= ACTION_ALIGNMENT_SKIP;
		act_			= 'alignment-deferred addition';
		alignDeferCount	= alignDeferCount+1;
		% Make a note so we don't try this next time
		% May be more than one in the model, which we need to note was
        % the cause of function 'nu' being rejected
		Aligned_out	= [Aligned_out ; nu*ones(numAligned,1)];
		Aligned_in	= [Aligned_in ; Used(findAligned)];
      end
    end
    %
    % Deletion: reinstate any previously deferred basis functions
    % resulting from this basis function
    % 
    if selectedAction==ACTION_DELETE
      findAligned	= find(Aligned_in==nu);
      numAligned	= length(findAligned);
      if numAligned>0
		reinstated					= Aligned_out(findAligned);
		Aligned_in(findAligned)		= [];
		Aligned_out(findAligned)	= [];
		%
		r_	= sprintf('%d ', reinstated);
		SB2_Diagnostic(OPTIONS,3,'Alignment reinstatement of %s\n', r_);
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% ACTION PHASE
  %%
  %% Implement above decision
  %%
  
  % We'll want to note if we've made a change which necessitates later
  % updating of the statistics
  % 
  UPDATE_REQUIRED	= false;
  %
  switch selectedAction
	%
   case ACTION_REESTIMATE,
    % 
    % Basis function 'nu' is already in the model, 
	% and we're re-estimating its corresponding alpha
	% 
	% - should correspond to Appendix A.3
    %
    oldAlpha	= Alpha(j);
    Alpha(j)	= newAlpha;
    s_j			= SIGMA(:,j);
    deltaInv	= 1/(newAlpha - oldAlpha);
    kappa		= 1/(SIGMA(j,j)+deltaInv);
    tmp			= kappa*s_j;
    SIGMANEW	= SIGMA - tmp * s_j';
    deltaMu		= -Mu(j)*tmp;
    Mu			= Mu + deltaMu;
    %
    if UpdateIteration
      S_in	= S_in + kappa*(BASIS_B_PHI * s_j).^2;
      Q_in	= Q_in - BASIS_B_PHI*(deltaMu);
    end
    updateCount	= updateCount+1;
    act_		= 're-estimation';
    %
    UPDATE_REQUIRED	= true;
    %
     
   case ACTION_ADD,
    % 
    % Basis function nu is not in the model, and we're adding it in
	% 
	% - should correspond to Appendix A.2
    %
    if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
      BASIS_Phi		= BASIS'*Phi;
      BASIS_PHI		= [BASIS_PHI BASIS_Phi];
      B_Phi			= beta*Phi;
      BASIS_B_Phi	= beta*BASIS_Phi;
    else
      B_Phi			= (Phi.*beta);
      BASIS_B_phi	= BASIS'*B_Phi;
    end
    tmp		= ((B_Phi'*PHI)*SIGMA)';
    %
    Alpha	= [Alpha ; newAlpha];
    PHI		= [PHI Phi]; 
    %
    s_ii		= 1/(newAlpha+S_in(nu));
    s_i			= -s_ii*tmp;
    TAU			= -s_i*tmp';
    SIGMANEW	= [SIGMA+TAU s_i ; s_i' s_ii];
    mu_i		= s_ii*Q_in(nu);
    deltaMu		= [-mu_i*tmp ; mu_i];
    Mu			= [Mu ; 0] + deltaMu;
	%
    if UpdateIteration
      mCi	= BASIS_B_Phi - BASIS_B_PHI*tmp;
      S_in	= S_in - s_ii * mCi.^2;
      Q_in	= Q_in - mu_i * mCi;
    end
    Used		= [Used; nu];
    addCount	= addCount+1;
    act_		= 'addition';
    %
    UPDATE_REQUIRED	= true;

   case ACTION_DELETE,
    % 
    % Basis function nu is in the model, but we're removing it
	% 
	% - should correspond to Appendix A.4
    %
    if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
      BASIS_PHI(:,j)	= [];
    end
    PHI(:,j)	= [];
    Alpha(j)	= [];
    %
    s_jj			= SIGMA(j,j);
    s_j				= SIGMA(:,j);
    tmp				= s_j/s_jj;
    SIGMANEW		= SIGMA - tmp*s_j';
    SIGMANEW(j,:)	= [];
    SIGMANEW(:,j)	= [];
    deltaMu			= - Mu(j)*tmp;
    mu_j			= Mu(j);
    Mu				= Mu +deltaMu;
    Mu(j)			= [];
	%
    if UpdateIteration
      jPm	= (BASIS_B_PHI * s_j);
      S_in	= S_in + jPm.^2 / s_jj;
      Q_in	= Q_in + jPm * mu_j / s_jj;
    end
    Used(j)		= [];
    deleteCount	= deleteCount+1;
    act_		= 'deletion';
    %
    UPDATE_REQUIRED	= true;
    %
  end					% of switch over actions
  M		= length(Used);
  %
  SB2_Diagnostic(OPTIONS, 3, 'ACTION: %s of %d (%g)\n', ...
				 act_, nu, deltaLogMarginal);
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% UPDATE STATISTICS
  
  % If we've performed a meaningful action,
  % update the relevant variables
  % 
  if UPDATE_REQUIRED
    %
    % S_in & Q_in values were calculated earlier
    % 
    % Here: update the S_out/Q_out values and relevance factors
    % 
    if UpdateIteration
	  %
	  % Previous "update" statisics calculated earlier are valid
	  % 
      S_out			= S_in;
      Q_out			= Q_in;
      tmp			= Alpha ./ (Alpha - S_in(Used));
      S_out(Used)	= tmp .* S_in(Used);
      Q_out(Used)	= tmp .* Q_in(Used);
      Factor		= (Q_out.*Q_out - S_out);
      SIGMA			= SIGMANEW;
      Gamma			= 1 - Alpha.*diag(SIGMA);
      if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian
		BASIS_B_PHI	= beta*BASIS_PHI;
      else
		BASIS_B_PHI	= ((PHI.* (beta*ones(1,M)))'*BASIS)';
      end
	else
	  %
	  % Compute all statistics in "full" form (non-Gaussian likelihoods)
	  %
      [SIGMA,Mu,S_in,Q_in,S_out,Q_out,Factor,newLogM,Gamma,...
       BASIS_B_PHI,beta] = ...
		SB2_FullStatistics(LIKELIHOOD,BASIS,PHI,Targets,Used,...
						   Alpha,beta,Mu,BASIS_PHI,BASIS_Targets,OPTIONS);
      deltaLogMarginal	= newLogM - logML;
    end
    %
    if UpdateIteration && deltaLogMarginal<0
      SB2_Diagnostic(OPTIONS,1,...
					 '** Alert **  DECREASE IN LIKELIHOOD !! (%g)\n',...
					 deltaLogMarginal);
    end
    %
    logML					= logML + deltaLogMarginal;
    count					= count + 1;
    logMarginalLog(count)	= logML;
  end
  
  % GAUSSIAN NOISE ESTIMATE
  % 
  % For Gaussian likelihood, re-estimate noise variance if:
  % 
  % - not fixed, AND
  % - an update is specified this cycle as normal, OR
  % - we're considering termination
  % 
  if LIKELIHOOD.InUse==LIKELIHOOD.Gaussian && ...
		~OPTIONS.fixedNoise && ...
		(selectedAction == ACTION_TERMINATE || ...
		 i<=CONTROLS.BetaUpdateStart || ...
		 rem(i,CONTROLS.BetaUpdateFrequency)==0)
	%
    betaZ1	= beta;
    y		= PHI * Mu;
    e		= (Targets-y);
    beta	= (N - sum(Gamma))/(e'*e);
    % Work-around zero-noise issue
    beta	= min([beta CONTROLS.BetaMaxFactor/var(Targets)]);
    %
    deltaLogBeta	= log(beta)-log(betaZ1);
	%
    if abs(deltaLogBeta) > CONTROLS.MinDeltaLogBeta
	  %
	  % Full re-computation of statistics required after beta update
	  % 
      [SIGMA,Mu,S_in,Q_in,S_out,Q_out,Factor,logML,Gamma,BASIS_B_PHI] = ...
		  SB2_FullStatistics(LIKELIHOOD,BASIS,PHI,Targets,Used,Alpha,beta,...
							 Mu, BASIS_PHI,BASIS_Targets,OPTIONS); 
      %
      count					= count + 1;
      logMarginalLog(count)	= logML;
      %
      if selectedAction==ACTION_TERMINATE
		%
		% We considered terminating above as no alpha update seemed
        % worthwhile. However, a beta update has made a non-trivial
        % increase in the likelihood, so we continue.
		% 
		selectedAction = ACTION_NOISE_ONLY;
		SB2_Diagnostic(OPTIONS,3,'Noise update (termination deferred)\n');
      end
    end
  end
  
  % CALLBACK
  %
  % Call callback function if specified
  % -	this can be useful for demos etc where it is desired to display
  %		graphical information at each iteration
  % 
  if OPTIONS.callback
    feval(OPTIONS.callbackFunc, i, selectedAction,...
		  logMarginalLog(1:count)/N, ...
		  Used, Mu./BasisScales(Used)', SIGMA,...
		  Alpha, beta, Gamma, PHI.*(ones(N,1)*BasisScales(Used)),...
		  OPTIONS.callbackData{:})
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % END OF CYCLE PROCESSING
  % 
  % Check if termination still specified, and output diagnostics
  % 
  if (selectedAction==ACTION_TERMINATE)
    %
    % If we're here, then no update to the model was considered worthwhile
    % 
    SB2_Diagnostic(OPTIONS,2,...
				   '** Stopping at iteration %d (Max_delta_ml=%g) **\n', ...
				   i, deltaLogMarginal);
    if LIKELIHOOD.InUse~=LIKELIHOOD.Gaussian
      SB2_Diagnostic(OPTIONS,2,'%4d> L = %.6f\t Gamma = %.2f (M = %d)\n',...
					 i, logML/N, sum(Gamma), M);
    else
      SB2_Diagnostic(OPTIONS,2,...
					 '%4d> L = %.6f\t Gamma = %.2f (M = %d)\t s=%.3f\n',...
					 i, logML/N, sum(Gamma), M, sqrt(1/beta));
    end
	% Exit the main loop
    break;
  end
  %
  % Check for "natural" termination
  %
  ITERATION_LIMIT	= (i==OPTIONS.iterations);
  TIME_LIMIT		= (etime(clock,t0)>OPTIONS.time);
  LAST_ITERATION	= ITERATION_LIMIT || TIME_LIMIT;
  %
  if (OPTIONS.monitor && ~rem(i,OPTIONS.monitor)) || LAST_ITERATION
	%
	% Output diagnostic progress info if desired
	% 
    if LIKELIHOOD.InUse~=LIKELIHOOD.Gaussian
      SB2_Diagnostic(OPTIONS,2,'%5d> L = %.6f\t Gamma = %.2f (M = %d)\n',...
		     i, logML/N, sum(Gamma), M);
    else
      SB2_Diagnostic(OPTIONS,2,...
		     '%5d> L = %.6f\t Gamma = %.2f (M = %d)\t s=%.3f \n',...
		     i, logML/N, sum(Gamma), M, sqrt(1/beta));
    end
  end
  
end % of MAIN LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% POST-PROCESSING
%%

%
% Warn if we exited the main loop without terminating automatically
% 
if (selectedAction~=ACTION_TERMINATE)
  %
  if ITERATION_LIMIT
	SB2_Diagnostic(OPTIONS,1,...
				   '** Iteration limit: algorithm did not converge\n');
  elseif TIME_LIMIT
	SB2_Diagnostic(OPTIONS,1,...
				   '** Time limit: algorithm did not converge\n');
  end
end
%
% Output action summary if desired
% 
if OPTIONS.diagnosticLevel>1
  % Stop timer
  t1	= etime(clock,t0);
  total	= addCount + deleteCount + updateCount;
  if CONTROLS.BasisAlignmentTest
    total	= total+alignDeferCount;
  end
  %
  SB2_Diagnostic(OPTIONS,2,'Action Summary\n');
  SB2_Diagnostic(OPTIONS,2,'==============\n');
  SB2_Diagnostic(OPTIONS,2,'Added\t\t%6d (%.0f%%)\n',...
		 addCount, 100*addCount/total);
  SB2_Diagnostic(OPTIONS,2,'Deleted\t\t%6d (%.0f%%)\n',...
		 deleteCount, 100*deleteCount/total);
  SB2_Diagnostic(OPTIONS,2,'Reestimated\t%6d (%.0f%%)\n',...
		 updateCount, 100*updateCount/total);
  %
  if CONTROLS.BasisAlignmentTest && alignDeferCount
    SB2_Diagnostic(OPTIONS,2,'--------------\n');
    SB2_Diagnostic(OPTIONS,2,'Deferred\t%6d (%.0f%%)\n',...
		   alignDeferCount, 100*alignDeferCount/total);
  end
  %
  SB2_Diagnostic(OPTIONS,2,'==============\n');
  SB2_Diagnostic(OPTIONS,2,'Total of %d likelihood updates\n', count);
  SB2_Diagnostic(OPTIONS,2,'Time to run: %s\n', SB2_FormatTime(t1));
end

% Terminate diagnostics
OPTIONS = SB2_Diagnostic(OPTIONS, 'end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% OUTPUT VARIABLES
%% 
if nargout>=1
  % We also choose to sort here - it can't hurt and may help later
  [PARAMETER.Relevant, index]	= sort(Used);
  % Not forgetting to correct for normalisation too
  PARAMETER.Value	= Mu(index) ./ BasisScales(Used(index))';
  if nargout>=2
	%
    HYPERPARAMETER.Alpha	= Alpha(index)./(BasisScales(Used(index))'.^2);
    HYPERPARAMETER.beta		= beta;
    if nargout>=3
	  %
      DIAGNOSTIC.Gamma		= Gamma(index);
      DIAGNOSTIC.Likelihood	= logMarginalLog(1:count);
      DIAGNOSTIC.iterations	= i;
      DIAGNOSTIC.S_Factor	= S_out;
      DIAGNOSTIC.Q_Factor	= Q_out;
      DIAGNOSTIC.M_full		= M_full;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
