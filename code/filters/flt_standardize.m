function [signal,state] = flt_standardize(varargin)
% Standardize a continuous EEG set causally.
% [Signal,State] = flt_standardize(Signal, State, WindowLength)
%
% Standardization ensures that per-channel data that normally can have any variance (e.g. due to
% varying conductivity, amplifier settings, etc.) is approximately normally distributed over data
% sets and time periods. This is usually necessary when learning and running models across sessions
% and subjects. It should not be applied before other artifact-rejection steps, as these steps
% usually take relative signal power into account. It is important to make the standardization
% window long enough that it does not factor out changes in signal power that one is interested in.
%
% Note that this function requires the data to be relatively free of artifacts to work well.
%
% In:
%   Signal       :   continuous data set to be filtered
%
%   WindowLength :   length of the window based on which normalization shall be performed, in
%                    seconds (default: 30)
%
%   State        :   previous filter state, as obtained by a previous execution of flt_iir on an
%                    immediately preceding data set (default: [])
%
% Out: 
%   Signal       :  standardized, continuous data set
%
%   State        :  state of the filter, after it got applied
%
% Examples:
%   % standardize the data in a moving window of default length (30s)
%   eeg = flt_standardize(eeg)
%
%   % standardize the data in a moving window of 60s length
%   eeg = flt_standardize(eeg,60)
%
%   % as previous, but passing all parameters by name
%   eeg = flt_standardize('Signal',eeg,'WindowLength',60)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17

% flt_standardize_version<1.01> -- for the cache

if ~exp_beginfun('filter') return; end

% follows IIR/FIR, as it should operate on a clean signal (rather than depend on HF noise, etc.)
declare_properties('name',{'Standardization','Standardize'}, 'cannot_follow','set_makepos', 'follows',{'flt_iir','flt_fir'}, 'independent_channels',true, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'window_len','WindowLength'}, 30, [], 'Length of normalization window. In seconds.'), ...
    arg({'multivariate_sphering','Sphere','sphere'}, false, [], 'Perform multivariate sphering. This removes correlations between channels but maintains. Assumes that the data is approximately zero-mean (i.e., first highpass filtered).'), ...
    arg({'stepsize','StepSize'}, 1/3, [], 'Step size between updates. The sphering matrix will be updated every this many samples. If this is below 1, it is assumed to be in seconds.','guru',true), ...
    arg({'usegpu','UseGPU'}, false, [], 'Use the GPU for sphering.'), ...
    arg({'lambda','CovarianceRegularization'}, 0.001, [], 'Covariance regularization. This is a regularization parameter for the covariance estimate used in sperhing.','guru',true), ...
    arg_nogui({'state','State'}));

warning off MATLAB:nearlySingularMatrix

% calc maximum amount of memory to use
if usegpu
    dev = gpuDevice(); maxmem = dev.FreeMemory/2^20;
else
    maxmem = hlp_memfree/(2^21);
end

% number of data points for our normalization window
N = round(window_len*signal.srate); %#ok<*NODEF>
if stepsize < 1
    stepsize = round(stepsize*signal.srate); end

% make up prior state if necessary
if isempty(state)
    for fld = utl_timeseries_fields(signal)
        field = fld{1};
        if size(signal.(field),2) > N
            % filter conditions & constant overall data offset (for better numerical accuracy; this is
            % unrelated to the running mean)
            state.(field) = struct('ord1',[],'ord2',[],'offset',sum(signal.(field) ,2)/size(signal.(field) ,2),'last_R',[]);
            % prepend a made up data sequence
            signal.(field) = [repmat(2*signal.(field) (:,1),1,N) - signal.(field)(:,(N+1):-1:2) signal.(field)];
        elseif ~isequal(signal.(field),1) && ~isempty(signal.(field))
            disp_once(['Not filtering the field .' field ': needs to be longer than the set data window length (for this data set ' num2str(window_len) ' seconds).']);
        end
    end
    prepended = true;
else
    prepended = false;
end

for fld = utl_timeseries_fields(signal)
    field = fld{1};
    if ~isempty(signal.(field)) && isfield(state,field)
        % get rid of NaN's and Inf's
        signal.(field)(~isfinite(signal.(field)(:))) = 0;
        
        if ~multivariate_sphering
            % get raw data X and running mean E[X]
            X = bsxfun(@minus,double(signal.(field)),state.(field).offset);
            [X_mean,state.(field).ord1] = moving_average(N,X,state.(field).ord1);
            % get squared data X^2 and running squared mean E[X^2]
            X2 = X.^2;
            [X2_mean,state.(field).ord2] = moving_average(N,X2,state.(field).ord2);
            % compute running std deviation sqrt(E[X^2] - E[X]^2)
            X_std = sqrt(abs(X2_mean - X_mean.^2));

            % compute standardized data
            signal.(field) = (X - X_mean) ./ X_std;
        else
            [C,S] = size(signal.(field));
            % split up the total sample range into k chunks that will fit in memory
            splits = ceil((C*C*S*8*8 + C*C*8*S/stepsize + C*S*8*2 + S*8*5) / (maxmem*1024*1024));
            if splits > 1
                fprintf('Now sphering data in %i blocks',splits); end            
            for i=1:splits
                range = 1+floor((i-1)*S/splits) : min(S,floor(i*S/splits));
                if ~isempty(range)
                    % get the data range
                    X = double(signal.(field)(:,range));
                    % move it to the GPU if applicable
                    if usegpu && length(range) > 1000
                        try X = gpuArray(X); catch,end; end
                    % compute running mean covariance (assuming a zero-mean signal)
                    [Xcov,state.(field).ord2] = moving_average(N,reshape(bsxfun(@times,reshape(X,1,C,[]),reshape(X,C,1,[])),C*C,[]),state.(field).ord2);
                    % extract the subset of time points at which we intend to update the sphering matrix
                    update_at = min(stepsize:stepsize:(size(Xcov,2)+stepsize-1),size(Xcov,2));
                    % if there is no previous R (from the end of the last chunk), we estimate it right at the first sample
                    if isempty(state.(field).last_R)
                        update_at = [1 update_at];
                        state.(field).last_R = eye(C);
                    end
                    Xcov = reshape(Xcov(:,update_at),C,C,[]);
                    if usegpu
                        Xcov = gather(Xcov); end
                    % do the reconstruction in intervals of length stepsize (or shorter if at the end of a chunk)
                    last_n = 0;
                    for j=1:length(update_at)
                        % update the sphering matrix
                        V = Xcov(:,:,j);
                        R = real((V*(1-lambda) + lambda*eye(C)*trace(V)/C)^(-1/2));
                        % apply the reconstruction to intermediate samples (using raised-cosine blending)
                        n = update_at(j);
                        subrange = range((last_n+1):n);
                        blend = (1-cos(pi*(1:(n-last_n))/(n-last_n)))/2;
                        signal.(field)(:,subrange) = bsxfun(@times,blend,R*signal.(field)(:,subrange)) + bsxfun(@times,1-blend,state.(field).last_R*signal.(field)(:,subrange));
                        [last_n,state.(field).last_R] = deal(n,R);
                    end
                end
                if splits > 1
                    fprintf('.'); end
            end
            if splits > 1
                fprintf('\n'); end            
        end
        
        % trim the prepended part if there was one
        if prepended
            signal.(field)(:,1:N) = []; end
    end
end

if usegpu
    state.ord2 = gather(state.ord2); end

exp_endfun;



function [X,Zf] = moving_average(N,X,Zi)
% Run a moving-average filter along the second dimension of the data.
% [X,Zf] = moving_average(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Channels x #Samples]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2012-01-10

if nargin <= 2 || isempty(Zi)
    Zi = zeros(size(X,1),N); end

% pre-pend initial state & get dimensions
Y = [Zi X]; M = size(Y,2);
% get alternating index vector (for additions & subtractions)
I = [1:M-N; 1+N:M];
% get sign vector (also alternating, and includes the scaling)
S = [-ones(1,M-N); ones(1,M-N)]/N;
% run moving average
X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% read out result
X = X(:,2:2:end);

if nargout > 1
    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end
