function targ = set_gettarget(signal)
% Generic function to extract the target values from a data set (epoched/continuous).
% Target = set_gettarget(Signal)
%
% Data sets may have associated "target values" in their meta-data (usually one per epoch and/or per
% event). These target values encode what outputs a cognitive monitoring system is supposed to
% output at the given epoch/marker time point/etc. These meta-data are the primary means in which
% the relationship between raw data and cognitive state is specified.
%
% If an epoched data set is passed, the .target field of the epochs will be used. If a continuous
% data set is passed which has events with a non-empty .target field, the target values of those
% events will be returned. Otherwise, if the signal has channels with a .type field that is set to
% 'target', the time course of these channels will be returned. Otherwise, the signal has no associated
% target values and this function will throw an error.
%
% In:
%   Signal  : EEGLAB data set, either epoched or continuous
%
% Out:
%   Target  : Array target values for the data set; either per epoch, per event, or per sample
%             this is of the form [N x D] where N is the number of target values and D is the dimensionality
%             of the target values (usually 1).
%
% Examples:
%   % obtain the target values for the given set
%   labels = set_gettarget(eeg)
%
% See also:
%   set_targetmarkers, set_makepos
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-07

if ~exist('signal','var')
    % trivial
    targ = [];
else
    % we got a bundle: extract target markers from the first stream
    if isfield(signal,'streams')
        signal = signal.streams{1}; end
    
    if all(isfield(signal,{'head','parts'})) && strcmp(char(signal.head),'set_partition') && length(signal.parts) >= 2
        % a computational shortcut applies if we're operating on partitioned data: this allows us to skip
        % the actual partitioning and instead look just at the target markers
        sourcetargets = set_gettarget(signal.parts{1});
        targ = sourcetargets(signal.parts{2},:);
    else
        
        % make sure that the data is evaluated
        signal = exp_eval_optimized(signal);
        
        if isfield(signal,'epoch') && ~isempty(signal.epoch) && isfield(signal.epoch,'target')
            % epoched data set with target field: get per-epoch targets
            targets = {signal.epoch.target};
            targetmask = ~cellfun('isempty',targets);
            if any(targetmask)
                targ = concat_targets(targets(targetmask));
                return;
            end
        end
        
        % continuous data set: check for events with non-empty target field
        if isfield(signal,'event') && isfield(signal.event,'target')
            targets = {signal.event.target};
            targetmask = ~cellfun('isempty',targets);
            if any(targetmask)
                if ~issorted([signal.event.latency])
                    warning('BCILAB:unsorted_events','The events in this data set are unsorted - this is likely an error in your processing pipeline.'); end
                targ = concat_targets(targets(targetmask));
                return;
            end
        end
        
        if isfield(signal,'chanlocs') && isfield(signal.chanlocs,'type') && isfield(signal,'data') && size(signal.data,3) == 1
            % continuous data set: get per-sample epoch targets
            targchans = strcmpi('target',{signal.chanlocs.type});
            if any(targchans)
                targ = signal.data(targchans,:)';
                return;
            end
        end
        
        % otherwise...
        error('set_gettarget did not find any target information in this data set. See help of set_gettarget and set_targetmarkers for how data sets can be annotated with target information.');
        
    end
end


function targets = concat_targets(targets)
if all(cellfun('size',targets,1) == 1)
    if length(unique(cellfun('size',targets,2))) > 1
        error('The target values in this set must have uniform dimensionality.'); end    
    targets = vertcat(targets{:});
elseif all(cellfun('size',targets,2) == 1)    
    if length(unique(cellfun('size',targets,1))) > 1
        error('The target values in this set must have uniform dimensionality.'); end    
    targets = horzcat(targets{:})';
else
    error('The target values in this set must be either all row vectors or all column vectors.');
end
