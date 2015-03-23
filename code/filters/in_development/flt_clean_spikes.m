function signal = flt_clean_spikes(varargin)
% Remove few-sample spikes from the data (by repeating the previous value).
% [Signal,State] = flt_clean_spikes(Signal, WindowLength, Cutoff, State)
%
% This method calculates the running median and running median absolute deviation of the signal in a
% sliding window, and flags each sample as a spike that is beyond a cutoff given in (robust) 
% standard deviations from the distribution in the sliding window. Eeach flagged sample is replaced 
% by the previous signal value.
%
% Note that, if spike removal is performed on a signal with large slow potential shifts (i.e., if
% applied prior to any high-pass filtering), sufficiently steep and strong trends in the signal can
% be flagged as spikes and replaced by previous signal values. Such data may prevent the use of
% aggressive thresholds.
%
% In:
%   Signal : Continuous EEGLAB dataset structs
%
%   Cutoff : Spike cutoff in standard deviations. This is in robust standard deviations of the
%            signal with respect to the moving window; any sample that lies outside this range is
%            considered a spike. If this is zero, the filter is equivalent to a median filter.
%            (default: 5)
%
%   WindowLength : Window length for the filter. In seconds. (default: 1)
%
%   State : optionally the filter state from a previous call
% 
% Out:
%   Signal : processed data
% 
%   State  :  state of the filter, after it got applied
%
% Examples:
%   % use the defaults
%   eeg = flt_clean_spikes(eeg)
%
%   % use a cutoff of 3 standard deviations
%   eeg = flt_clean_spikes(eeg, 3)
%
%   % use a window length of 0.5 seconds and a cutoff of 7 standard deviations
%   eeg = flt_clean_spikes('Signal',eeg, 'WindowLength',0.5, 'Cutoff',7)
%
%   % run a basic median filter
%   eeg = flt_clean_spikes('Signal',eeg,'Cutoff',0)
%
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-11-06

% flt_clean_spikes_version<0.8> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name','SpikeRemoval', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'stddev_cutoff','Cutoff','StandardDevCutoff','cutoff'}, 5, [0 0 15 1000], 'Spike cutoff in standard deviations. This is in robust standard deviations of the signal with respect to the moving window; any sample that lies outside this range is considered a spike. If this is zero, the filter is equivalent to a median filter.'), ...
    arg({'window_len','WindowLength'}, 1, [0 10000], 'Window length for the filter. In seconds.'), ...
    arg_nogui({'state','State'}));

% number of data points for our normalization window
N = round(window_len*signal.srate); %#ok<*NODEF>


% for each time series field...
for f = utl_timeseries_fields(signal)
    field = f{1};
    if ~isempty(signal.(field))
        % flip dimensions so that we can filter along the 1st dimension
        [X,dims] = spatialize_transpose(double(signal.(field)));

        % get rid of NaN's and Inf's
        X(~isfinite(X(:))) = 0;        

        % initialize filter state if necessary
        if ~isfield(state, field)
            % 1st and 2nd order filter conditions, last sample value
            state.(field) = struct('ord1',[],'ord2',[],'last',X(1,:));
            % we prepend the signal with a mirror section of itself, to minimize start-up transients
            % (and if the signal is too short, we repeat it as much as we need)
            X = [repmat(2*X(1,:),N,1) - X(1+mod(((N+1):-1:2)-1,size(X,1)),:); X]; %#ok<AGROW>
            prepended = true;
        else
            prepended = false;
        end
        
        if length(X) > 100000
            fprintf('flt_clean_spikes: applying running median; this may take a while...'); end

        % get running median E[X]            
        [X_med,state.(field).ord1] = running_median(N,X,state.(field).ord1);
        
        if stddev_cutoff > 0
            % get absolute difference between running median and X
            X_diff = abs(X - X_med);
            % get running median of that (running median absolute deviation)
            [X_mad,state.(field).ord2] = running_median(N,X_diff,state.(field).ord2);
            % flag spikes
            spike_mask = X_diff > X_mad*1.4826*stddev_cutoff;
            % for each channel...
            for c=1:size(X,2)
                % interpolate each spike by holding the previous value
                for k=find(spike_mask(:,c))'
                    if k==1
                        X(k,c) = state.(field).last(c);
                    else
                        X(k,c) = X(k-1,c);
                    end
                end
            end
        else
            % basic median filter if cutoff is 0
            X = X_med;
        end
        
        if length(X) > 100000
            fprintf('done.\n'); end

        % cut off the data segment that we previously prepended
        if prepended
            X(1:N,:) = []; end

        % unflip dimensions and write the result back
        signal.(field) = unspatialize_transpose(X,dims);
    end
end

exp_endfun;


function [X,Zf] = running_median(N,X,Zi)
% Apply a running median filter along the first dimension of X.
% [X,Zf] = running_median(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Samples x #Channels]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2014-11-06

% ensure that the window length is odd
if ~mod(N,2)
    N = N+1; end
N_half = floor(N/2);

% handle initial and final state
if nargin <= 2 || isempty(Zi)
    Zi = zeros(N-1,size(X,2)); end
X = [Zi; X];
if nargout > 1
    Zf = X(end-N+2:end,:); end

% apply running median for each channel
for c=1:size(X,2)
    X(N:end,c) = fastmedfilt1d(X(N_half+1:end-N_half,c),N,X(1:N_half,c),X(end-N_half+1:end)); end

X = X(1:end-N+1,:);