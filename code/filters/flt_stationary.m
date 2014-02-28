function signal = flt_stationary(varargin)
% Restrict an epoched signal to an n-dimensional stationary subspace
% Signal = flt_stationary(Signal, StationaryDim)
%
% This is an implementation of the Analytical Stationary Subspace Analysis algorithm [1]. Typically,
% it is applied to continuous data with a block length that is large enough to average over changes
% in the experimental conditions (if presented in a randomized bock design), e.g. 30 seconds.
%
% The data typically needs to be high-pass filtered prior to use, since otherwise drifts will
% dominate the results. The filter may also be applied to epoched data, in which case the epochs are 
% taken as the blocks.
%
% In:
%   Signal        : continuous data set to be filtered
%
%   StationaryDim : assumed number of stationary sources; if a negative number is given, it is taken
%                   as the assumed number of non-stationary sources; may also be given as a ratio,
%                   in which case it is understood as a fraction of the number of channels (default: -0.15)
%
%   BlockLength   : length of the signal blocks across which non-stationarity should be assessed,
%                   in seconds; may also be an [Nx2] array of intervals that should be taken as blocks
%                   (default: 30)
%
%   Operation     : Operation to perform; can be one of the following:
%                   * 'keep_stationary': project the signal onto the stationary components
%                   * 'keep_nonstationary': project the signal onto the non-stationary components
%                   * 'separate': order the signal by increasing non-stationarity of its components
%                   * 'backproject_stationary': back-project the stationary components onto the channels
%                   * 'backproject_nonstationary': back-project the non-stationary components onto the channels
%
% Out:
%   Signal        :  filtered, EEGLAB data set
%
% Examples:
%   % use default settings to remove the 10 least stationary components of the signal
%   eeg = flt_stationary(eeg);
%
%   % retain a 20-dimensional maximally stationary subspace, as measured across blocks of 15 seconds length
%   eeg = flt_stationary(eeg,20,15)
%
%   % as before, but passing arguments by name
%   eeg = flt_stationary('Signal',eeg,'StationaryDim',20,'BlockLength',15)
%
%   % as before, but this time retain all except for a 20-dimensional maxially non-stationary subpace 
%   eeg = flt_stationary('Signal',eeg,'StationaryDim',-20,'BlockLength',15)
%
%   % as before, but now back-project the retained subspace onto the channels (giving a rank-deficit signal)
%   eeg = flt_stationary('Signal',eeg,'StationaryDim',-20,'BlockLength',15,'Operation','backproject_stationary')
%
% References:
%  [1] S. Hara, Y. Kawahara, P. von Buenau, "Stationary Subspace Analysis as a Generalized Eigenvalue Problem",
%      Lecture Notes in Computer Science, 2010, Vol. 6443/2010, 422-429.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-07-28

% flt_stationary_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% we prefer to run this after spectral filters (e.g. FIR), and before the epoch extraction
declare_properties('name','StationarySubspace', 'precedes',{'set_makepos','flt_fir'}, 'follows','flt_iir',  'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'numstat','StationaryDim'}, -0.15, [], 'Assumed # of stationary sources. If a negative number is given, it is taken as the assumed number of non-stationary sources.'), ...
    arg({'blocklen','BlockLength'}, 30, [0 Inf], 'Block length. Length of the signal blocks across which non-stationarity should be assessed, in seconds. May also be an [Nx2] array of intervals that should be taken as blocks.'), ...
    arg({'opmode','Operation'}, 'backproject_stationary',{'keep_stationary','keep_nonstationary','separate','backproject_stationary','backproject_nonstationary'}, 'Operation to perform.'), ...
    arg_norep('decomposition',unassigned),...
    arg_norep('new_chanlocs',unassigned));

if ~exist('decomposition','var')
    % we need to compute the subspace decomposition first        
    
    % create the blocks X{i}
    if ndims(signal.data) == 3
        % signal is already epoched; ignore blocklen
        X = squeeze(mat2cell(signal.data,signal.nbchan,signal.pnts,ones(1,signal.trials)));
    elseif isscalar(blocklen)
        % generate blocklen-sized blocks
        l = blocklen*signal.srate;
        X = mat2cell(signal.data,signal.nbchan,[l*ones(1,floor(signal.pnts/l)), mod(signal.pnts,l)]);
        X = X(cellfun('size',X,2) == l);
    else
        % extract intervals for each row in blocklen...
        for i=size(blocklen,1):-1:1
            X{i} = signal.data(:,round(blocklen(i,1)*signal.srate):round(blocklen(i,2)*signal.srate)); end
    end
    
    N = length(X);    
    % compute mean and covariance for each block
    for i=N:-1:1
        mu{i} = mean(X{i}'); %#ok<UDIM>
        sig{i} = cov(X{i}');
    end
    
    % and compute joint mean and covariance
    Mu = mean(vertcat(mu{:}));
    Sig = mean(cat(3,sig{:}),3);
    invSig = inv(Sig);
    
    % compute the matrix S (Eq. 9)
    S = zeros(size(Sig));
    for i=1:N
        S = S + mu{i}*mu{i}' + (1/2) * sig{i} * invSig * sig{i}; end %#ok<MINV>
    S = S/N - Mu*Mu' - 1/2*Sig;
    
    % solve the generalized eigenvalue problem and sort results
    [phi,lambdas] = eig(S,Sig);
    [lambdas,idx] = sort(diag(lambdas),'ascend'); %#ok<ASGLU>
    phi = phi(:,idx);
    
    % split into stationary and non-stationary subspaces
    if abs(numstat) < 1
        numstat = round(numstat*signal.nbchan); end
    if numstat < 0
        numstat = signal.nbchan+numstat; end        
    stationary = phi(:,1:numstat)';
    nonstationary = phi(:,numstat+1:end)';
    
    switch opmode
        case 'keep_stationary'
            decomposition = stationary;
        case 'keep_nonstationary'
            decomposition = nonstationary;
        case 'separate'
            decomposition = phi';
        case 'backproject_stationary'
            decomposition = stationary' * stationary;
        case 'backproject_nonstationary'
            decomposition = nonstationary' * nonstationary;
        otherwise
            error('Unsupported operation requested.');
    end
    
    if ~any(strcmp(opmode,{'backproject_stationary','backproject_nonstationary'}))        
        new_chanlocs = struct('labels',cellfun(@(x)sprintf('StationaryComponent%.0f',x),num2cell(1:size(decomposition,1),1),'UniformOutput',false)); 
    else
        new_chanlocs = [];
    end
end

% project data
[C,S,T] = size(signal.data); %#ok<*NODEF>
signal.data = reshape(decomposition*reshape(signal.data,C,[]),[],S,T);
signal.nbchan = size(signal.data,1);

% rewrite chanlocs if necessary
if ~isempty(new_chanlocs)
    signal.chanlocs = new_chanlocs; end

exp_endfun('append_online',{'decomposition',decomposition,'new_chanlocs',new_chanlocs});