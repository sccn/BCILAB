function signal = flt_clean_channels(varargin)
% Remove channels with abnormal data from a continuous data set.
% Signal = flt_clean_channels(Signal,MinCorrelation,IgnoredQuantile,WindowLength,MaxBrokenTime,Rereferenced)
%
% This is an automated artifact rejection function which ensures that the data contains no channels
% that record only noise for extended periods of time. If channels with control signals are
% contained in the data these are usually also removed. There are two threshold criteria: one is a
% minimum required correlation between a channel and a surrogate of it calculated from its neighbors
% using spline interpolation (calculated in a manner that is robust to bad channels in the
% neighborhood) and the other is a maximum tolerated noise level in standard deviations relative to
% the remaining channels (also robust).
%
% In:
%   Signal          : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%                     with a 0.5Hz - 2.0Hz transition band).
%
%   CorrelationThreshold : Correlation threshold. If a channel is correlated at less than this value
%                          to its reconstruction from other channels, it is considered abnormal in
%                          the given time window. Note that this method can only be used when
%                          channel locations are available. (default: 0.8)
%
%   LineNoiseThreshold : If a channel has more line noise relative to its signal than this value, in
%                        standard deviations from the channel population mean, it is considered abnormal.
%                        (default: 4)
%
%
%   The following are "detail" parameters that usually do not have to be tuned. If you can't get
%   the function to do what you want, you might consider adapting these to your data.
%   
%   NumSamples : Number of RANSAC samples. This is the number of samples to generate in the random
%                sampling consensus process. (default: 50)
%
%   SubsetSize : Subset size. This is the size of the channel subsets to use, as a fraction of the
%                total number of channels. (default: 0.25)
%
%   WindowLength    : Length of the windows (in seconds) for which correlation is computed; ideally
%                     short enough to reasonably capture periods of global artifacts or intermittent 
%                     sensor dropouts, but not shorter (for statistical reasons). (default: 5)
% 
%   MaxBrokenTime : Maximum time (either in seconds or as fraction of the recording) during which a 
%                   retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                   (very lax). (default: 0.4)
%
%	ProtectChannels : list of channel names (cell array) that should be protected from removal. 
%                     (default: {})
%
%
% The following arguments are deprecated but retained for backwards compatibility:
%
%   Rereferenced    : Whether the measures should be computed on re-referenced data. This can improve 
%                     performance in environments with extreme EM noise, but will decrease robustness 
%                     against individual channels with extreme excursions. (default: false)
%
%   LineNoiseAware : Whether the operation should be performed in a line-noise aware manner. If enabled,
%                    the correlation measure will not be affected by the presence or absence of line 
%                    noise. (default: true).
%
%   MinCorrelation  : Minimum correlation between a channel and any other channel (in a short period 
%                     of time) below which the channel is considered abnormal for that time period.
%                     Reasonable range: 0.4 (very lax) to 0.6 (quite aggressive). (default: 0.5). 
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
%                                2014-05-12

% flt_clean_channels_version<0.9.7> -- for the cache

if ~exp_beginfun('filter') return; end;

declare_properties('name','ChannelCleaning', 'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'corr_threshold','CorrelationThreshold'}, 0.8, [0 0.3 0.95 1], 'Correlation threshold. If a channel is correlated at less than this value to its reconstruction from other channels, it is considered abnormal in the given time window. Note that this method can only be used when channel locations are available.'), ...
    arg({'noise_threshold','LineNoiseThreshold'},4,[],'Line-noise threshold. If a channel has more line noise relative to its signal than this value, in standard deviations from the channel population mean, it is considered abnormal.'), ...
    arg({'window_len','WindowLength'}, 5, [0 0.25 5 Inf], 'Window length to compute correlations. Length of the windows (in seconds) for which correlation is computed; ideally short enough to reasonably capture periods of global artifacts (which are ignored), but not shorter (for statistica reasons).'), ...
    arg({'max_broken_time','MaxBrokenTime','ignored_time','MaxIgnoredTime'}, 0.4, [0 Inf], 'Maximum duration/fraction of broken data to tolerate. Maximum time (either in seconds or as fraction of the recording) during which a retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6 (very lax).'), ...
    arg({'num_samples','NumSamples'}, 50, [], 'Number of RANSAC samples. This is the number of samples to generate in the random sampling consensus process.','guru',true), ...
    arg({'subset_size','SubsetSize'}, 0.25, [], 'Subset size. This is the size of the channel subsets to use, as number of channels or a fraction of the total number of channels.','guru',true), ...
    arg({'protect_channels','ProtectChannels'},[],[],'Channels to protect from removal. This protects the channels with the given names from being removed.','type','cellstr','shape','row'), ...
    arg({'keep_unlocalized_channels','KeepUnlocalizedChannels'},false,[],'Keep unlocalized channels. Whether to keep channels which have no localiztion information and can therefore not be checked based on location information.'), ...
    arg({'use_gpu','UseGPU'}, false, [], 'Whether to run on the GPU. Makes sense for offline processing if you have a GTX Titan or better.'), ...
    arg_deprecated({'linenoise_aware','LineNoiseAware'},true,[],'Line-noise aware processing. Whether the operation should be performed in a line-noise aware manner. If enabled, the correlation measure will not be affected by the presence or absence of line noise.','guru',true), ...
    arg_deprecated({'rereferenced','Rereferenced'},false,[],'Run calculations on re-referenced data. This can improve performance in environments with extreme EM noise, but will decrease robustness against individual channels with extreme excursions.'), ...
    arg_deprecated({'min_corr','MinimumCorrelation'}, 0.5, [0 1], 'Minimum correlation between channels. If the measure falls below this threshold in some time window, the window is considered abnormal.'), ...
    arg_deprecated({'ignored_quantile','IgnoredQuantile'}, 0.1, [0 1], 'Quantile of highest correlations ignored. Upper quantile of the correlation values that may be arbitrarily high without affecting the outcome - avoids problems with shorted channels.'), ...
    arg_norep('removed_channel_mask',unassigned)); 

% flag channels
if ~exist('removed_channel_mask','var')
    subset_size = round(subset_size*size(signal.data,1)); 

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

    if linenoise_aware && signal.srate > 100
        % remove signal content above 50Hz
        B = design_fir(100,[2*[0 45 50]/signal.srate 1],[1 1 0 0]);
        for c=signal.nbchan:-1:1
            X(:,c) = filtfilt_fast(B,1,signal.data(c,:)'); end
        % determine z-scored level of EM noise-to-signal ratio for each channel
        noisiness = mad(signal.data'-X,1)./mad(X,1);
        znoise = (noisiness - median(noisiness)) ./ (mad(noisiness,1)*1.4826);        
        % trim channels based on that
        noise_mask = znoise > noise_threshold;
    else
        X = signal.data';
        noise_mask = false(C,1);
    end

    % optionally subtract common reference from data
    if rereferenced
        X = bsxfun(@minus,X,mean(X,2)); end
    
    if isfield(signal.chanlocs,'X') && isfield(signal.chanlocs,'Y') && isfield(signal.chanlocs,'Z') && all([length([signal.chanlocs.X]),length([signal.chanlocs.Y]),length([signal.chanlocs.Z])] > length(signal.chanlocs)*0.5)
        fprintf('Scanning for bad channels...');
        
        % get the matrix of all channel locations [3xN]
        [x,y,z] = deal({signal.chanlocs.X},{signal.chanlocs.Y},{signal.chanlocs.Z});
        usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
        locs = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
        X = X(:,usable_channels);
        
        P = hlp_diskcache('filterdesign',@calc_projector,locs,num_samples,subset_size);
        corrs = zeros(length(usable_channels),W);
        
        if use_gpu
            P = gpuArray(P);
            X = gpuArray(X);
            corrs = gpuArray(corrs);
        end
        
        % calculate each channel's correlation to its RANSAC reconstruction for each window
        for o=1:W
            XX = X(offsets(o)+wnd,:);
            YY = sort(reshape(XX*P,length(wnd),length(usable_channels),num_samples),3);
            YY = YY(:,:,round(end/2));
            corrs(:,o) = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
        end
        
        if use_gpu
            corrs = gather(corrs); end
        flagged = corrs < corr_threshold;
        
        % mark all channels for removal which have more flagged samples than the maximum number of
        % ignored samples
        removed_channel_mask = quickif(keep_unlocalized_channels,false(C,1),true(C,1));
        removed_channel_mask(usable_channels) = sum(flagged,2)*window_len > max_broken_time;
    else
        fprintf('No channel locations were available; falling back to determining bad channels based on their mutual correlation...');    
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
    end
    
    % also incorporate the line noise criterion
    removed_channel_mask = removed_channel_mask(:) | noise_mask(:);
    
    fprintf(' removing %i channels...\n',nnz(removed_channel_mask));
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


% calculate a bag of reconstruction matrices from random channel subsets
function P = calc_projector(locs,num_samples,subset_size)
stream = RandStream('mt19937ar','Seed',435656);
rand_samples = {};
for k=num_samples:-1:1
    tmp = zeros(size(locs,2));
    subset = randsample(1:size(locs,2),subset_size,stream);
    tmp(subset,:) = real(sphericalSplineInterpolate(locs(:,subset),locs))';
    rand_samples{k} = tmp;
end
P = horzcat(rand_samples{:});


function Y = randsample(X,num,stream)
Y = [];
while length(Y)<num
    pick = round(1 + (length(X)-1).*rand(stream));
    Y(end+1) = X(pick);
    X(pick) = [];
end

function Y = mad(X,usemedians)
if usemedians
    Y = median(abs(bsxfun(@minus,X,median(X))));
else
    Y = mean(abs(bsxfun(@minus,X,mean(X))));
end
