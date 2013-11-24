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
%   EventTypes : Cell array of event valuelist (strings) that encode a prediction target value 
%                (i.e. a condition of interest).
%
%                The order in which marker valuelist are supplied determines the corresponding associated
%                target value of each. For example, if this is {'xxx' 'dfg' 'yyy'}, all events with 
%                type 'yyy' will be assigned a target value of 3, while those with type 'dfg' will 
%                receive target value of 2, etc. (i.e., the target value is the index of the 
%                event type, as it appears in the cell array).
%
%                To assign the same target value to events of multiple valuelist, group them into a
%                nested cell array, e.g., {'xxx' {'dfg' 'blah'} 'yyy'} -- here, the assignment
%                from event type to target value is as follows: xxx -> 1, dfg -> 2, blah -> 2, yyy -> 3.
%
%   EventMap : Alternative specification of the mapping from event type to target value; given 
%              as a cell array of {'type1',value1,'type2',value2,'type3',value3} where the respective
%              type denotates an event type, and the value is the target value that will be associated
%              with each event of the respective type.
%
%   EventField : Field of the EEGLAB event structure which contains the events to be annotated, provided as
%                a string. If no event field is given, the 'type' field is used by default.
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
    arg({'eventmap','EventMap','eventtypes','EventTypes'}, '', [], 'Marker type to target value map. This is either a cell array of ''type'', (numeric) value pairs, or a cell array of strings/cellstrings where the strings in the k''th cell element are mapped onto value k. Can also be a single cell containing the special string ''actualvalues'', which automatically maps each type to its numeric value.','shape','row','type','expression'), ...
    arg({'epoch_bounds','EpochBounds'}, [0 0], [], 'Assumed epoch boundaries. How much data around each marker must be present for it to be considered a potential target marker.','shape','row'), ...
    arg({'eventfield','EventField'}, 'type', [], 'Field of event structure containing the events.'), ...
    arg({'prune_nontarget','PruneNontarget'}, false, [], 'Prune non-target events. Whether to prune non-target events from the data.'), ...
    arg({'avoid_boundaries','AvoidBoundaries'}, true, [], 'Avoid boundary events. Whether potential target markers whose epoch bounds would overlap markers with type ''boundary'' should be excluded.'));

if ~isempty(signal.event)
    types = {signal.event.(eventfield)};
    matchidx = zeros(size(types));
    if isequal(eventmap,{'actualvalues'})
        matchidx = ~cellfun('isempty',types);
    else
        % translate eventmap to typelist/valuelist
        chars = cellfun('isclass',eventmap,'char');
        cells = cellfun('isclass',eventmap,'cell');
        numbers = cellfun('isclass',eventmap,'double') | cellfun('isclass',eventmap,'single');
        if (all(chars(1:2:end)) && all(numbers(2:2:end)) && length(eventmap)>1 && mod(length(eventmap),2)==0)
            % {'type',value,'type',value,...} format
            typelist = eventmap(1:2:end);
            valuelist = eventmap(2:2:end);
        elseif all(chars|cells)
            % {'type',{'type','type'},'type','type', ...} format
            typelist = {}; valuelist = {};
            for i=1:length(eventmap)
                if ~iscell(eventmap{i})
                    typelist{end+1} = eventmap{i}; %#ok<*AGROW>
                    valuelist(end+1) = {i};
                else
                    typelist(end+(1:length(eventmap{i}))) = eventmap{i};
                    valuelist(end+(1:length(eventmap{i}))) = {i};
                end
            end
        else
            error('BCILAB:set_targetmarkers:unrecognized_syntax','Unrecognized syntax for EventMap parameter: %s',hlp_tostring(eventmap));
        end
        % perform matching
        for t=1:length(typelist)
            if any(typelist{t}=='?' | typelist{t}=='*')
                matchidx(~cellfun('isempty',regexp(types,['^',strrep(strrep(typelist{t},'?','.'),'*','.{0,}'),'$']))) = t;
            else
                matchidx(strcmp(types,typelist{t})) = t;
            end
        end
    end
    
    % get candidate target event indices
    candidates = find(matchidx ~= 0);
    if ~isempty(candidates)
        % prune candidates that exceed the data set boundaries
        latencies = round([signal.event(candidates).latency]);
        eporange = round(epoch_bounds(1)*signal.srate) : round(epoch_bounds(2)*signal.srate);
        candidates(latencies+eporange(1)<1 | latencies+eporange(end)>signal.pnts | latencies<1 | latencies>signal.pnts) = [];

        % optionally prune candidates whose epochs cross boundary events
        if avoid_boundaries && isempty(signal.epoch) && ~isempty(candidates)
            boundaries = strcmp({signal.event.type},'boundary');
            if any(boundaries)
                % generate a sparse mask of boundary events
                boundlats = min(signal.pnts,max(1,round([signal.event(boundaries).latency])));
                boundmask = sparse(ones(1,length(boundlats)),boundlats,1,1,signal.pnts);
                % prune events that intersect the boundary mask
                latencies = round([signal.event(candidates).latency]);
                if ~isempty(latencies)
                    candidates(any(boundmask(bsxfun(@plus,eporange',latencies)))) = []; end
            end
        end
    end
    
    % assign the target values for the retained events
    if ~isempty(candidates)
        if isequal(eventmap,{'actualvalues'})
            [signal.event(candidates).target] = celldeal(cellfun(@str2num,types(candidates),'UniformOutput',false));
        else
            [signal.event(candidates).target] = celldeal(valuelist(matchidx(candidates))); 
        end
    elseif ~isfield(signal.event,'target')
        signal.event(1).target = [];
    end
    
    % optionally prune non-target events
    if prune_nontarget
        signal.event(~matchidx) = []; end
end

% update misc fields
signal.etc.epoch_bounds = epoch_bounds;

exp_endfun;
