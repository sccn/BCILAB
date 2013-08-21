function signal = set_targetmarkers(varargin)
% Associate a subset of markers with a target value for prediction.
% Signal = set_targetmarkers(Signal,EventTypes,EventMap,EventField,PruneNonTarget)
%
% This function sets up an association between events in a data set (with a specific type) on one
% hand and values that a BCI is supposed to predict (called "target values"). It is the primary way
% by which the desired output of a BCI for some data is defined. Note: this function is typically
% called by bci_train (according to its TargetMarkers parameter).
%
% After events have been annotated with this function, a BCI approach calibrated on this data set
% may extract epochs of data around these "target" events and learn a mapping that predicts the
% respective target value whenever it encounters data effectively similar to the segments/epochs
% that had that target value associated with them in the calibration data. This way, events that
% denote a particular user condition in the data can be assigned a concrete value that a BCI should
% output whenever the user is in that condition again.
%
% The function effectively adds a .target field to the events in the data set, and sets it to desired
% value for those events that are selected in either EventTypes or EventMap. The remaining events 
% will have an empty target field.
%
% In:
%   Signal     : continuous EEGLAB data set which shall be annotated with target markers
%
%   EventTypes : Cell array of event types (strings) that encode a prediction target value 
%                (i.e. a condition of interest).
%
%                The order in which marker types are supplied determines the corresponding associated
%                target value of each. For example, if this is {'xxx' 'dfg' 'yyy'}, all events with 
%                type 'yyy' will be assigned a target value of 3, while those with type 'dfg' will 
%                receive target value of 2, etc. (i.e., the target value is the index of the 
%                event type, as it appears in the cell array).
%
%                To assign the same target value to events of multiple types, group them into a
%                nested cell array, e.g., {'xxx' {'dfg' 'blah'} 'yyy'} -- here, the assignment
%                from event type to target value is as follows: xxx -> 1, dfg -> 2, blah -> 2, yyy -> 3.
%
%   EventMap : Alternative specification of the mapping from event type to target value; given 
%              as a cell array of {'type1',value1,'type2',value2,'type3',value3} where the respective
%              type denotates an event type, and the value is the target value that will be associated
%              with each event of the respective type.
%
%   EventField : Field of the EEGLAB event structure which contains the events to be annotated, provided as
%              a string. If no event field is given, the 'type' field is used by default.
%
%   PruneNontarget: Prune non-target events. Whether to prune non-target events from the data.
%
% Out:
%   Signal : data set with "target" annotations added to the events of interest (.target field added)
%            (continuous or epoched)
%
% Examples:
%   % annotate all events of type 'TX' with target value 1 and all events of type 'TY' with target value 2
%   eeg = set_targetmarkers(eeg,{'TX','TY'})
%   
%   % as before, but assign target value 1 to 'TY' events and 2 to 'TX' events
%   eeg = set_targetmarkers(eeg,{'TY','TX'})
%
%   % as before, but assign target value 1 to both 'TX' and 'TY' events, assign target value 2
%   % to 'TZ' events, and 3 to 'A' and 'B' events.
%   eeg = set_targetmarkers(eeg,{{'TY','TX'},'TZ',{'A','B'}})
%
%   % as before, but using named arguemnts
%   eeg = set_targetmarkers('Signal',eeg, 'EventTypes',{{'TY','TX'},'TZ',{'A','B'}})
%
%   % as before, but express it in the alternative EventMap specification
%   eeg = set_targetmarkers('Signal',eeg, 'EventMap',{'A',3, 'B',3, 'TX',1, 'TY',1, 'TZ',2})
%
%   % use some fairly arbitrary target values this time
%   eeg = set_targetmarkers('Signal',eeg, 'EventMap',{'A',0.75, 'B',-2, 'TX',1000,'TY',1001})
%
% See also:
%   set_makepos, bci_train, set_gettarget
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-27

% set_targetmarkers_version<1.1> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('name',{'TargetMarkers'},'independent_channels',true,'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'eventtypes','EventTypes'}, {}, [], 'Marker types to be annotated. Cell array of strings. The order in which marker types are supplied determines the target value that is assigned to each respective event (the assigned value is the index within this list).','type','cellstr','shape','row'), ...
    arg({'eventmap','EventMap'}, {}, [], 'Marker type / target value map. Cell array of ''type'',value pairs, else a single cell containing the string ''actualvalues'', which automatically maps each type to its numeric value. This is an alternative specification of how event types map onto target values.','shape','row'), ...
    arg({'epoch_bounds','EpochBounds'}, [], [], 'Assumed epoch boundaries. How much data around each marker must be present for it to be considered a potential target marker.','shape','row'), ...
    arg({'eventfield','EventField'}, 'type', [], 'Field of event structure containing the events.','type','char','shape','row'), ...
    arg({'prune_nontarget','PruneNontarget'}, false, [], 'Prune non-target events. Whether to prune non-target events from the data.'));

typelist = {};
targetlist = {};
if ~isempty(eventtypes)
    % create a flattened eventtype list and a target value list
    if any(cellfun(@iscell,eventtypes))
        for i=1:length(eventtypes)
            if ~iscell(eventtypes{i})
                targetlist(end+1) = {i};
                typelist{end+1} = eventtypes{i};
            else
                targetlist(end+1:end+length(eventtypes{i})) = {i};
                typelist(end+1:end+length(eventtypes{i})) = eventtypes{i};
            end
        end           
    else
        typelist = eventtypes;
        targetlist = num2cell(1:length(eventtypes));
    end
end

if ~isempty(eventmap)
    % append to the typelist and targetlist
    
    % use the actual numeric values of the events as their targets
    if length(eventmap) == 1 && strcmp(eventmap{1}, 'actualvalues')
        if isempty(signal.event)
            warning('BCILAB:set_targetmarkers:no_events',['The given signal has no events: ' hlp_tostring(signal)]);  
        else
            allstring = {signal.event.(eventfield)};
            emptyvals = cellfun(@isempty, allstring);
            allstring(emptyvals) = {'NaN'};
            allnumeric = cellfun(@str2double, allstring);
            numerictypes = ~isnan(allnumeric); 
            typelist = [typelist allstring(numerictypes)];
            targetlist = [targetlist num2cell(allnumeric(numerictypes))];
        end
    else % event map is type, value pairings
        for k=2:2:length(eventmap)
            if ~isnumeric(eventmap{k})
                error('EventMap is a cell array of the form ''type'',value,''type'',value, ... where all values must be numeric.'); end
        end
        typelist = [typelist eventmap(1:2:end)];
        targetlist = [targetlist eventmap(2:2:end)];
    end
end


if isempty(signal.event)
    warning('BCILAB:set_targetmarkers:no_events',['The given signal has no events: ' hlp_tostring(signal)]);    
else
    if ~isfield(signal.event,'target')
        signal.event(1).target = []; end
    
    if isempty(signal.epoch)
        if ~issorted([signal.event.latency])
            disp('set_targetmarkers: Events in this data set are unsorted; sorting them now.');
            [sorted,idx] = sort([signal.event.latency]); %#ok<ASGLU>
            signal.event = signal.event(idx);
        end
    end
    
    % make a mapping from event idx onto the entry # in the targetlist (or zero)
    alltypes = {signal.event.(eventfield)};
    matchidx = zeros(size(alltypes));
    for t=1:length(typelist)
        if any(typelist{t}=='?' | typelist{t}=='*')
            matchidx(~cellfun(@isempty,regexp(alltypes,['^',strrep(strrep(typelist{t},'?','.'),'*','.{0,}'),'$']))) = t;
        else
            matchidx(cellfun(@(x)isequal(x,typelist{t}),alltypes)) = t;
        end
    end
    
    % optionally prune non-target events
    if prune_nontarget
        signal.event(~matchidx) = [];
        matchidx(~matchidx) = [];
    end
    
    % get a mask of retained events
    retain = find(matchidx ~= 0);
    if isempty(retain)
        disp('Note: None of the events in the data set corresponded to any of the target markers.'); end
    
    % assumes epoch bounds are okay if data is already epoched
    if ~isempty(epoch_bounds) && isempty(signal.epoch)
        % epoch bounds were given: further restrict the set of retained events
        % generate epoch index range, in samples
        eporange = round(epoch_bounds(1)*signal.srate) : round(epoch_bounds(2)*signal.srate);
        
        if ~isempty(eporange)
            % prune events that exceed the data set boundaries
            lats = round([signal.event(retain).latency]);
            retain(lats+eporange(1)<1 | lats+eporange(end)>signal.pnts) = [];
            
            % generate a sparse mask of boundary events
            boundlats = min(signal.pnts,max(1,round([signal.event(strcmp({signal.event.type},'boundary')).latency])));
            if ~isempty(boundlats)
                boundmask = sparse(ones(1,length(boundlats)),boundlats,1,1,signal.pnts);
                
                % prune events that intersect the boundary mask
                lats = round([signal.event(retain).latency]);
                if ~isempty(lats)
                    retain(any(boundmask(bsxfun(@plus,eporange',lats)))) = []; end
            end
        end
    end
    
    % finally assign the target values for the retained events
    if ~isempty(retain)
        [signal.event(retain).target] = celldeal(targetlist(matchidx(retain)));
    else
        disp('Note: no target markers were selected.');
    end
    
    signal.etc.epoch_bounds = epoch_bounds;
end


exp_endfun;
