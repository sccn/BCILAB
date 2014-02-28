function signal = flt_clean_channels(varargin)
% Remove channels with abnormal data from a continuous data set.
% Signal = flt_clean_channels(Signal,MinCorrelation,IgnoredQuantile,WindowLength,MaxBrokenTime,Rereferenced)
%
% This is an automated artifact rejection function which ensures that the data contains no channels
% that record only noise for extended periods of time. If channels with control signals are
% contained in the data these are usually also removed. The criterion is based on correlation: if a
% channel is decorrelated from all others (pairwise correlation < a given threshold), excluding a
% given fraction of most correlated channels -- and if this holds on for a sufficiently long fraction 
% of the data set -- then the channel is removed.
%
% In:
%   Signal          : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%                     with a 0.5Hz - 2.0Hz transition band).
%
%   MinCorrelation  : Minimum correlation between a channel and any other channel (in a short period 
%                     of time) below which the channel is considered abnormal for that time period.
%                     Reasonable range: 0.4 (very lax) to 0.6 (quite aggressive). (default: 0.5). 
%                     
%
%   The following are "detail" parameters that usually do not have to be tuned. If you can't get
%   the function to do what you want, you might consider adapting these to your data.
%   
%   IgnoredQuantile : Fraction of channels that need to have at least the given MinCorrelation value
%                     w.r.t. the channel under consideration. This allows to deal with channels or
%                     small groups of channels that measure the same noise source, e.g. if they are
%                     shorted. If many channels can be disconnected during an experiment and you
%                     have strong noise in the room, you might increase this fraction, but consider
%                     that this a) requires you to decrease the MinCorrelation appropriately and b)
%                     can make the correlation measure more brittle. Reasonable range: 0.05 (rather
%                     lax) to 0.2 (very tolerant re disconnected/shorted channels).The default is
%                     0.1.
%
%   WindowLength    : Length of the windows (in seconds) for which correlation is computed; ideally
%                     short enough to reasonably capture periods of global artifacts (which are
%                     ignored), but not shorter (for statistical reasons). (default: 2)
% 
%   MaxBrokenTime : Maximum time (either in seconds or as fraction of the recording) during which a 
%                   retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                   (very lax). (default: 0.4)
%
%   Rereferenced    : Whether the measures should be computed on re-referenced data. This can improve 
%                     performance in environments with extreme EM noise, but will decrease robustness 
%                     against individual channels with extreme excursions. (default: false)
%
%   LineNoiseAware : Whether the operation should be performed in a line-noise aware manner. If enabled,
%                    the correlation measure will not be affected by the presence or absence of line 
%                    noise. (default: true).
%
%	ProtectChannels : list of channel names (cell array) that should be protected from removal. 
%                     (default: {})
%
% Out:
%   Signal : data set with bad channels removed
%
% Examples:
%   % use with defaults
%   eeg = flt_clean_channels(eeg);
%
%   % override the MinimumCorrelation and the IgnoredQuantile defaults
%   eeg = flt_clean_channels(eeg,0.7,0.15);
%
%   % override the MinimumCorrelation and the MaxIgnoredTime, using name-value pairs
%   eeg = flt_clean_channels('Signal',eeg,'MinimumCorrelation',0.7, 'MaxBrokenTime',0.15);
%
%   % override the MinimumCorrelation and the MaxIgnoredTime, using name-value pairs 
%   % in their short forms
%   eeg = flt_clean_channels('signal',eeg,'min_corr',0.7, 'max_broken_time',0.15);
%
% See also:
%   flt_clean_windows, flt_clean_settings
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-07-06

% flt_clean_channels_version<0.8.2> -- for the cache

if ~exp_beginfun('filter') return; end;

declare_properties('name','ChannelCleaning', 'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'min_corr','MinimumCorrelation'}, 0.5, [0 1], 'Minimum correlation between channels. If the measure falls below this threshold in some time window, the window is considered abnormal.'), ...
    arg({'ignored_quantile','IgnoredQuantile'}, 0.1, [0 1], 'Quantile of highest correlations ignored. Upper quantile of the correlation values that may be arbitrarily high without affecting the outcome - avoids problems with shorted channels.'), ...
    arg({'window_len','WindowLength'}, 2, [], 'Window length to compute correlations. Length of the windows (in seconds) for which correlation is computed; ideally short enough to reasonably capture periods of global artifacts (which are ignored), but not shorter (for statistica reasons).'), ...
    arg({'max_broken_time','MaxBrokenTime','ignored_time','MaxIgnoredTime'}, 0.4, [], 'Maximum duration/fraction of broken data to tolerate. Maximum time (either in seconds or as fraction of the recording) during which a retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6 (very lax).'), ...
    arg_deprecated({'rereferenced','Rereferenced'},false,[],'Run calculations on re-referenced data. This can improve performance in environments with extreme EM noise, but will decrease robustness against individual channels with extreme excursions.'), ...
    arg({'linenoise_aware','LineNoiseAware'},true,[],'Line-noise aware processing. Whether the operation should be performed in a line-noise aware manner. If enabled, the correlation measure will not be affected by the presence or absence of line noise.','guru',true), ...
    arg({'protect_channels','ProtectChannels'},{},[],'Channels to protect from removal. This protects the channels with the given names from being removed.','type','cellstr','shape','row'), ...
    arg_norep('removed_channel_mask',unassigned)); 

% flag channels
if ~exist('removed_channel_mask','var')
    if max_broken_time > 0 && max_broken_time < 1  %#ok<*NODEF>
        max_broken_time = size(signal.data,2)*max_broken_time;
    else
        max_broken_time = signal.srate*max_broken_time;
    end
    
    [C,S] = size(signal.data);
    window_len = window_len*signal.srate;
    wnd = 0:window_len-1;
    offsets = round(1:window_len:S-window_len);
    W = length(offsets);    
    retained = 1:(C-ceil(C*ignored_quantile));
        
    % optionally ignore both 50 and 60 Hz spectral components...
    if linenoise_aware && signal.srate > 110
        if signal.srate <= 130
            B = design_fir(500,[2*[0 45 50 55]/signal.srate 1],[1 1 0 1 1]);
        else
            B = design_fir(500,[2*[0 45 50 55 60 65]/signal.srate 1],[1 1 0 1 0 1 1]);
        end
        for c=signal.nbchan:-1:1
            X(:,c) = filtfilt_fast(B,1,signal.data(c,:)'); end
    else
        X = signal.data';
    end

    % optionally subtract common reference from data
    if rereferenced
        X = bsxfun(@minus,X,mean(X,2)); end
    
    % for each window, flag channels with too low correlation to any other channel (outside the
    % ignored quantile)
    flagged = zeros(C,W);
    for o=1:W
        sortcc = sort(abs(corrcoef(X(offsets(o)+wnd,:))));
        flagged(:,o) = all(sortcc(retained,:) < min_corr);
    end
    
    % mark all channels for removal which have more flagged samples than the maximum number of
    % ignored samples
    removed_channel_mask = sum(flagged,2)*window_len > max_broken_time;
    fprintf('Removing %i channels...\n',nnz(removed_channel_mask));
    
    % remove the channels in the protect list
    if ~isempty(protect_channels)
        removed_channel_mask(set_chanid(signal,protect_channels)) = true; end    
end

% annotate the data with what was removed (for visualization)
if ~isfield(signal.etc,'clean_channel_mask')
    signal.etc.clean_channel_mask = true(1,signal.nbchan); end
signal.etc.clean_channel_mask(signal.etc.clean_channel_mask) = ~removed_channel_mask;

% execute
if any(removed_channel_mask)
    signal.data = signal.data(~removed_channel_mask,:,:,:,:,:,:,:);
    signal.chanlocs = signal.chanlocs(~removed_channel_mask);
    signal.nbchan = size(signal.data,1);
end

exp_endfun('append_online',{'removed_channel_mask',removed_channel_mask});
