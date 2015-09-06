function signal = flt_clean_flatlines(varargin)
% Remove channels with abnormal data from a continuous data set.
% Signal = flt_clean_flatlines(Signal,MaxFlatlineDuration,MaxAllowedJitter)
%
% This is an automated artifact rejection function which ensures that 
% the data contains no flat-lined channels.
%
% In:
%   Signal : continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%            with a 0.5Hz - 2.0Hz transition band)
%
%   MaxFlatlineDuration : Maximum tolerated flatline duration. In seconds. If a channel has a longer
%                         flatline than this, it will be considered abnormal. Default: 5
%
%   MaxAllowedJitter : Maximum tolerated jitter during flatlines. As a multiple of epsilon.
%                      Default: 20
%
% Out:
%   Signal : data set with flat channels removed
%
% Examples:
%   % use with defaults
%   eeg = flt_clean_flatlines(eeg);
%
% See also:
%   flt_clean_windows, flt_clean_peaks
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-07-06

% flt_clean_flatlines_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end;

declare_properties('name','FlatlineCleaning', 'follows','flt_selchans', 'precedes','flt_laplace', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'max_flatline_duration','MaxFlatlineDuration'}, 5, [0 Inf], 'Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal.'), ...
    arg({'max_allowed_jitter','MaxAllowedJitter'}, 20, [0 Inf], 'Maximum tolerated jitter during flatlines. As a multiple of epsilon.'), ...
    arg_deprecated({'min_stddev','MinimumStdDev'}, 0.0001, [0 Inf], 'This parameter has no function any more.'), ...
    arg_norep('removed_channels',unassigned)); 

% flag channels
if ~exist('removed_channels','var')
    removed_channels = [];
    for c=1:signal.nbchan
        zero_intervals = reshape(find(diff([false abs(diff(signal.data(c,:)))<(max_allowed_jitter*eps) false])),2,[])';
        if max(zero_intervals(:,2) - zero_intervals(:,1)) > max_flatline_duration*signal.srate
            removed_channels(end+1) = c; end
    end
end

% execute
if ~isempty(removed_channels)
    retain_channels = true(1,size(signal.data,1)); 
    retain_channels(removed_channels) = false;
    signal.data = signal.data(retain_channels,:,:,:,:,:,:,:);
    signal.chanlocs = signal.chanlocs(retain_channels);
    signal.nbchan = size(signal.data,1);
end

exp_endfun('append_online',{'removed_channels',removed_channels});
