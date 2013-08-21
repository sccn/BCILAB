function signal = set_makepos(varargin)
% Create epochs relative to the target markers in a data set.
% Signal = set_makepos(Signal,TimeWindow)
%
% This function turns a continuous data set with "target" markers (aka target events), i.e. events
% that have a .target field which contains a non-empty "target value" into an epoched data set,
% creating one epoch relative to each target marker, and copying that markers's .target value into the
% respective epoch's .target field. Target markers are created via the function set_targetmarkers.
%
% In:
%   Signal      :   continuous data set from which epochs shall be extracted
%
%   TimeWindow  :   time window relative to the events (forming the extracted epoch), in seconds
%
% Out:
%   Signal  :  epoched data set, with target variable assigned; changes:
%              - data is reformatted into a 3d array
%              - xmin & xmax indicate the epoch bounds in seconds
%              - events are rewritten to be relative to epochs (urevents unaffected)
%              - times contains the time points inside the epoch
%
% Notes:
%   In the regular BCILAB processing pipeline, target markers are defined on the continuous 
%   data set either automatically by the bci_ functions, or manually via set_targetmarkers.
%
% Examples:
%   % extract epochs from a continuous data set, beginning 0.5s before the position of each marker
%   % and ending 1.5s after the position of the corresponding marker
%   eeg = set_makepos(eeg,[-0.5 1.5])
%
%   % as before, but first add target markers for event types 'X' and 'Y'
%   eeg = set_makepos(set_targetmarkers(eeg,{'X','Y'}),[-0.5 1.5])
%
% See also:
%   set_targetmarkers, bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_makepos_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name',{'EpochExtraction','epoch'}, 'independent_channels',true,'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'epobounds','TimeWindow'}, [], [], 'Epoch time window relative to the target markers. In seconds.','type','expression','shape','row'));

if ~onl_isonline
    evtbackup = signal.event;
    had_icaact = ~isempty(signal.icaact) && ~isscalar(signal.icaact) && ~isempty(signal.icaweights) && ~isempty(signal.icasphere) && ~isempty(signal.icachansind);
    if isempty(signal.event)
        error('This data set has no events.'); end
    
    % make sure that the events are sorted (otherwise pop_epoch may behave erratically)
    if ~issorted([signal.event.latency])
        disp('set_makepos: The events in this data set are unsorted; sorting them now.');
        [sorted,idx] = sort([signal.event.latency]); %#ok<ASGLU>
        signal.event = signal.event(idx);
    end
    
    % identify the target events
    if ~isfield(signal.event,'target')
        % If you are getting this error, then the events of your signal are lacking a piece of
        % meta-data, namely the field signal.event.target. The regular processing pipeline (bci_train
        % or the model calibration GUI) automatically add this field, so either you are calling this
        % function directly on "raw" data (in which case you should first run set_targetmarkers on it),
        % or this field got lost during one of the intermediate processing stages (in which case
        % there is an error in one of the filter plugins or the basic EEGLAB functions).
        error('The events of this signal must have a .target field (which is assigned using set_targetmarkers).'); 
    end
    targmask = ~cellfun('isempty',{signal.event.target});
    
    if ~any(targmask)
        % Note: if you are getting this error it means that no event in the data was tagged as a
        % "target" event (these are the ones that the BCI is supposed to deal with). The tagging
        % happens automatically during bci_train / "calibrate model" based on the TargetMarkers that
        % you specified, so the most likely reason for this error is that no event type in the data
        % matched any of those that you gave as "target markers". Please check whether you specified
        % them in accordance with what's in your data. Another, less likely, reason is that you
        % accidentally "lost" them during some manual data curation or during a pre-processing step
        % (unlikely, unless you are using a lot of ad hoc plugins).
        error('This data set contains none of the target events (i.e. events with non-empty .target field). Did you specify the correct target markers?'); 
    end
    
    % temporarily back up the .type field into .oldtype and replace the type of the target events by '__target__'
    [signal.event.oldtype] = signal.event.type;
    [signal.event(targmask).type] = deal('__target__');

    % extract epochs & do sanity checks
    [signal, evtindices] = pop_epoch(signal, {'__target__'}, epobounds);
    if length(evtindices) ~= size(signal.data,3) || ((length(evtindices) ~= length(signal.epoch)) && ~isempty(signal.epoch))
        % Note: If you are getting this error, the event information returned by pop_epoch is not 
        % what is expected by BCILAB. This should not happen, but pop_epoch is a very complex function;
        % therefore, if you get this message, you probably spotted a bug deep in the plumbing of 
        % EEGLAB -- please report it.
        error('BCILAB:set_makepos:inconsistent_epochinfo','The data returned by pop_epoch is inconsistent; unable to determine the per-epoch target value.');
    end
    
    if length(evtindices) ~= nnz(targmask)
        % Note: If you are getting this error, there are some events in your data set that have been
        % labeled as "target" events (by set_targetmarkers) but which are so close to a boundary
        % event or the edge of the dataset, that no epochs could be extracted for them. This may
        % either be caused by some intermediate processing stage that introduced new boundary events
        % into the data (e.g. artifact rejection) without properly pruning the invalidated target
        % markers, or, somewhat less likely, because set_targetmarkers was not given a large enough
        % "safety range" around the events (e.g., either you explicitly passed a range that is
        % smaller than the size of your actual epochs or bci_train failed to figure out the size of
        % your epochs by itself -- which would be a bug). Or, least likely, it may have been caused
        % by a bug already in set_targetmarkers (should not happen since the code is fairly
        % simple, but if you manage to track this down, please report it).
        error('BCILAB:set_makepos:inconsistent_issue','Failed to extract an epoch for each target event.');
    end
        
    % copy back the .oldtype field to the .type field, if some events are left
    if ~isempty(signal.event)
        [signal.event.type] = signal.event.oldtype;
        signal.event = rmfield(signal.event,'oldtype');
    end
    
    % make sure that the .epoch field has the right length
    % (it can be empty because sometimes no events are retained)
    if length(signal.epoch) ~= nnz(targmask)
        if isempty(signal.epoch)
            signal.epoch = struct('target',cell(1,nnz(targmask)));
        else
            signal.epoch(nnz(targmask)).target = [];
        end
    end
    
    % copy back old types to the epoch field
    if ~isempty(signal.epoch)
        try
            [signal.epoch.eventtype] = signal.epoch.eventoldtype;
            signal.epoch = rmfield(signal.epoch,'eventoldtype');
        catch
            % (... which might perhaps be malformed or non-exitant)
        end
    end
    
    % assign target values and latencies to the epoch
    [signal.epoch.target] = evtbackup(targmask).target;
    [signal.epoch.latency] = evtbackup(targmask).latency;
    
    % reconstruct the ICA activation if required
    if had_icaact        
        signal.icaact = reshape((signal.icaweights*signal.icasphere)*signal.data(signal.icachansind,:), [], signal.pnts, signal.trials); end
else
    % in the online case we don't extract epochs but instead behave as an epoch that slides over the data
    
    % fix up xmin and xmax for the current epoch
    signal.xmin = epobounds(1);
    signal.xmax = epobounds(1) + (size(signal.data,2)-1)/signal.srate;    
    % add an epoch field; this contains at least the latency (of the time-locking event, in samples)
    signal.epoch = struct('latency', {signal.smax - signal.xmax*signal.srate});        
end

exp_endfun;
