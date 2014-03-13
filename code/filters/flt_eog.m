function [signal,state] = flt_eog(varargin)
% Remove artifacts from EEG using artifact reference channels.
% [Signal,State] = flt_eog(Signal, ReferenceChannels, ForgetFactor, KernelLength, RemoveReferenceChannels)
%
% This is an online filter that operates on continuous data, and removes artifacts using a
% regression technique, if artifact channels (e.g., EOG or EMG) are present (using recursive least
% squares) [1]. Note that noise in the artifact signals may be transferred onto the EEG channels.
%
% In:
%   Signal :  continuous data set to be filtered
%
%   ReferenceChannels : list of artifact reference channel indices or cell array of channel names
%                       (default: [] = try to auto-detect)
%
%   ForgetFactor : forgetting factor of the adaptive filter; amounts to a choice of the 
%                  effective memory length (default: 0.9995)
%
%   KernelLength : length/order of the temporal FIR filter kernel (default: 3)
%
%   RemoveReferenceChannels : whether to remove the reference channels after processing (default: false)
%
%   State : previous filter state, as obtained by a previous execution of flt_eog on an
%           immediately preceding data set (default: [])
%
% Out:
%   Signal : filtered, continuous EEGLAB data set
%
%   State : state of the filter, can be used to continue on a subsequent portion of the data
%
%
% Examples:
%   % using the defaults
%   eeg = flt_eog(eeg)
%
%   % manually supply EOG channels
%   eeg = flt_eog(eeg,{'veog','heog'});
%
%   % pass a specific forgetting factor (here: by name)
%   eeg = flt_eog('Signal',eeg,'ForgetFactor',0.9995)
%
% References:
%  [1] P. He, G.F. Wilson, C. Russel, "Removal of ocular artifacts from electro-encephalogram by adaptive filtering"
%      Med. Biol. Eng. Comput. 42 pp. 407-412, 2004
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17

% flt_eog_version<1.03> -- for the cache

if ~exp_beginfun('filter') return; end

% makes no sense on epoched data (and should precede the major spatial filters)
declare_properties('name','EOGRemoval', 'precedes','flt_project', 'cannot_follow',{'set_makepos','flt_ica'}, 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'eogchans','ReferenceChannels','EOGChannels'}, [], [], 'Reference channels. These are the channels that carry the artifacts (e.g., EOG).','type','expression','shape','row'), ...
    arg({'ffact','ForgetFactor'}, 0.9995, [0.99 1], 'Forgetting factor. Determines the memory length of the adaptive filter.','guru',true), ...
    arg({'kernellen','KernelLength'}, 3, uint32([1 1 20 1000]), 'Kernel Length. The length/order of the temporal FIR filter kernel.'), ...
    arg({'removeeog','RemoveReferenceChannels','RemoveEOG'}, true, [], 'Remove reference channels. Remove the reference channels after processing.'), ...
    arg_nogui({'state','State'}));

if size(signal.data,3) > 1
    error('flt_eog is supposed to be applied to continuous (non-epoched) data.'); end

% initialize the state, if necessary
if isempty(state)
    % figure out what the EOG channels are
    if isempty(eogchans)
        eogchans = find(strcmp({signal.chanlocs.type},'EOG'));
    else
        eogchans = set_chanid(signal,eogchans);
    end
    if isempty(eogchans)
        disp_once('flt_eog: No artifact/EOG channels found in the data.'); end
    state.eog = eogchans;                          % eog channel indices
    state.eeg = setdiff(1:signal.nbchan,eogchans); % eeg channel indices
    state.neog = length(state.eog);                % number of eog channel indices
    
    % initialize RLS filter state
    state.hist = zeros(state.neog,kernellen);     % hist is the block of the M last eog samples in matrix form
    state.R_n = eye(state.neog * kernellen) / 0.01; % R(n-1)^-1 is the inverse matrix
    state.H_n = zeros(state.neog*kernellen,length(state.eeg));  % H(n-1) is the EOG filter kernel
end

% apply filter
[signal.data,state.hist,state.H_n,state.R_n] = compute(signal.data,state.hist,state.H_n,state.R_n,state.eeg,state.eog,ffact);

if removeeog
    removed_channel_mask = true(1,size(signal.data,1));
    removed_channel_mask(state.eeg) = false;
    % annotate the data with what was removed (for visualization)
    if ~isfield(signal.etc,'clean_channel_mask')
        signal.etc.clean_channel_mask = true(1,signal.nbchan); end
    signal.etc.clean_channel_mask(signal.etc.clean_channel_mask) = ~removed_channel_mask;
    % perform removal
    if any(removed_channel_mask) 
        signal.data = signal.data(~removed_channel_mask,:,:,:,:,:,:,:);
        signal.chanlocs = signal.chanlocs(~removed_channel_mask);
        signal.nbchan = size(signal.data,1);
    end
end

exp_endfun;


function [X,hist,H_n,R_n] = compute(X,hist,H_n,R_n,eeg,eog,ffact)
% for each sample...
for n=1:size(X,2)
    % update the EOG history by feeding in a new sample
    hist = [hist(:,2:end) X(eog,n)];
    % vectorize the EOG history into r(n)        % Eq. 23
    tmp = hist';
    r_n = tmp(:);
    
    % calculate K(n)                             % Eq. 25
    K_n = R_n * r_n / (ffact + r_n' * R_n * r_n);
    % update R(n)                                % Eq. 24
    R_n = ffact^-1 * R_n - ffact^-1 * K_n * r_n' * R_n;
    
    % get the current EEG samples s(n)
    s_n = X(eeg,n);    
    % calculate e(n/n-1)                         % Eq. 27
    e_nn = s_n - (r_n' * H_n)';    
    % update H(n)                                % Eq. 26
    H_n = H_n + K_n * e_nn';
    % calculate e(n), new cleaned EEG signal     % Eq. 29
    e_n = s_n - (r_n' * H_n)';
    % write back into the signal
    X(eeg,n) = e_n;
end
