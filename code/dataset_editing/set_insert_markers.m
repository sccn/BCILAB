function newsignal = set_insert_markers(varargin)
% Inject runs of markers into specific segments of a continuous data set.
% Signal = set_insert_markers(Signal, Options...)
%
% Almost all real-time inference in BCIs is done on the basis of epochs, and epochs are most
% conveniently created relative to certain ("time-locking") events/markers in the data. This
% function allows to cover periods of continuous data with events, at regular or random intervals,
% so that epochs covering these ranges can subsequently be extracted. What periods shall be
% populated with events can be flexibly specified.
%
% In:
%   Signal      : continuous data set
%
%   SegmentSpec : segment specification. cell array of one of the following forms (lats in seconds):
%                 (default: {0 Inf})
%                 note: the ordering of time values w.r.t. event-type values and cell-array values in the subsequent 
%                       specifications is arbitrary, whereas the ordering of time values w.r.t. each other is relevant 
%                       (first time value shall be lower than second time value); the ordering of eventtypes w.r.t. each 
%                       other is also relevant.
%                 * {absolute_time absolute_time}:
%                   segment specified using two time points
%                 * {event_type relative_time relative_time} / {relative_time  event_type relative_time} / 
%                   {relative_time relative_time event_type}:
%                   here, the segment is relative to some event (of a given type), in between the time interval given by the 
%                   first and second time value
%                 * {event_type relative_time relative_time event_type} / {relative_time event_type event_type relative_time} / ...
%                   here, the segment is in between two immediately successive events of the given types (any intervals with other 
%                   events in between the specified ones are not considered for injection), constrained by the relative lats 
%                   for each event
%                 * {event_type relative_time {ignore_type1,ignore_type2,...} relative_time event_type}
%                   as above, except that intermediate events of type ignore_type1/ignore_type2/etc. are ignored
%                   (if the in-between cell array is empty, any other events are ignored)
%
%   Limits : optional time limits (in seconds) to constrain event placement (default: [-Inf Inf])
%
%   Event : the inserted event type string, or alternatively a template event struct 
%           (default: "mytype_i", when injected relative to an event of type "mytype" or
%           "mytype1_mytype2", when injected in between two events of type "mytype1" and "mytype2")
%           * if a string is specified, all event fields besides the 'type', 'latency' and 'duration' fields 
%             contained in the data will be left empty 
%           * if a struct is specified, the 'latency' field will be substituted appropriately
%
%   Count : number of events inserted within an interval; see Counting for the counting scheme (default: 1)
%
%   Counting : what count means for any given segmentspec, either 'perinterval' or 'persecond' (default: 'perinterval')
%
%   Placement : how the injected events should be placed, either 'random' or 'equidistant' (default: 'equidistant')
%
%   Repeatable : whether the randomization procedure shall give repeatable results (default: 1); different numbers (aside from 0)
%                give different repeatable runs, i.e. the value determines the randseed
%
%   MinLength : segments that are shorter than this (in seconds) are ignored. (default: 0)
%
%   MaxLength : segments that are longer than this (in seconds) are ignored. (default: Inf)
%
% Out:
%   Signal  : continuous data set with new events injected
%
% Notes:
%   The only parameter that may be specified by position (instead of as name-value pair) is the first one.
%
% Examples:
%   % place 20 events of type 'X' within the interval 1000s to 2000s into the given data set (regular placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{1000 2000},'Count',20,'Event','X')
%
%   % as before, but use an event struct
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{1000 2000},'Count',20,'Event',struct('type','X'))
%
%   % place 20 events within the interval 1000s to 2000s into the given data set (random pleacement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{1000 2000},'Count',20,'Placement','random','Event','X')
%
%   % place 3 events per second within the interval 1000s to 2000s into the given data set (regular placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{1000 2000}, 'Count',3, 'Counting','persecond','Event','X')
% 
%   % place on average 3 events per second within the interval 1000s to 2000s into the given data set (random placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{1000 2000}, 'Count',3, 'Counting','persecond','Placement','random','Event','X')
%
%   % place 20 events of type 'X' within each interval within 2s to 10s following each occurrence of the event 'A'
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{'A',2,10},'Count',20,'Event','X')
%
%   % place on average 5 events per second (typed 'X') within each interval within -5s to 10s around each occurrence of the event 'A' (random placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{'A',-5,10},'Counting','persecond','Count',5,'Event','X','Placement','random')
%
%   % same as before, equivalent SegmentSpec formatting
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A',10},'Counting','persecond','Count',5,'Event','X','Placement','random')
%
%   % place 10 events (typed 'X') between each successive occurrence of event 'A' followed by event 'B' (with no other event in between),
%   % and begin the interval 5s after event 'A' and end it 3s before event 'B' (regular placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{'A',5,-3,'B'},'Count',10,'Event','X')
%
%   % place 10 events (typed 'X') between each successive occurrence of event 'A' followed by event 'B' (with no other event in between),
%   % and begin the interval 5s *before* event 'A' and end it right on event 'B' (regular placement)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A',0,'B'},'Count',10,'Event','X')
%
%   % as before, but insert 3 events per second
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A',0,'B'},'Count',3,'Counting','persecond','Event','X')
%
%   % as before, but also consider those intervals where other events of type 'p' and/or 'q' occur between the 'A' and the 'B'
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {'p','q'}, 0,'B'},'Count',3,'Counting','persecond','Event','X')
% 
%   % as before, but also consider intervals where any other event occurs between the 'A' and the 'B' (except for 'B' obviously)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {}, 0,'B'},'Count',3,'Counting','persecond','Event','X')
%
%   % as before, but discard segments that would be longer than 10 seconds
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {}, 0,'B'},'Count',3,'Counting','persecond','Event','X','MaxLength',10)
%
%   % as before, but use random placement
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {}, 0,'B'},'Count',3,'Counting','persecond','Event','X','Placement','random')
%
%   % as before, but use a random rand seed to obtain different placing at every call
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {}, 0,'B'},'Count',3,'Counting','persecond','Event','X','Placement','random','Repeatable',0)
%
%   % as before, but use a fixed specific rand seed to obtain a specific (but repeatable placing)
%   eeg = set_insert_markers(eeg, 'SegmentSpec',{-5,'A', {}, 0,'B'},'Count',3,'Counting','persecond','Event','X','Placement','random','Repeatable',10)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-24

% set_insert_markers_version<1.0> -- for the cache

global tracking;

if ~exp_beginfun('editing') return; end

declare_properties('name','MarkerInsertion', 'independent_channels',true, 'cannot_follow','set_makepos', 'independent_channels',true,'independent_trials',true);

opts = arg_define([0 1],varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_subswitch({'segmentspec','SegmentSpec','segment'},'absoluterange', ...
        {'absoluterange',{ ...
            arg({'lo','BeginOffset'},0,[],'Lower bound of insertion interval. Events will be inserted beginning from this time point, in seconds.'), ...
            arg({'hi','EndOffset'},Inf,[],'Upper bound of insertion interval. Events will be inserted up to this time point, in seconds. If this is negative, it counts from the end of the recording.')}, ...
         'relativerange',{ ...
            arg({'event','EventType'},'event1',[],'Reference event type. New events will be inserted in a range around each event of this type.'), ...
            arg({'lo','BeginOffset'},-0.5,[],'Lower bound relative to event. This is the lower boundary of insertion intervals relative to the reference events. In seconds.'), ...
            arg({'hi','EndOffset'},1,[],'Upper bound relative to event. This is the upper boundary of insertion intervals relative to the reference events. In seconds.')}, ...
         'spannedrange',{ ...
            arg({'openevent','OpenEvent'},'event1',[],'Type of opening reference event. New events will be inserted between each pair of successive events with types OpenEvent and CloseEvent.'), ...
            arg({'closeevent','CloseEvent'},'event2',[],'Type of closing reference event. New events will be inserted between each pair of successive events with types OpenEvent and CloseEvent.'), ...
            arg({'lo','OpenOffset'},0.5,[],'Offset relative to opening event. This is an offset relative to the position of the opening reference event, which shifts the beginning of a spanned insertion interval. In seconds.'), ...
            arg({'hi','CloseOffset'},-0.5,[],'Offset relative to closing event. This is an offset relative to the position of the closing reference event, which shifts the beginning of a spanned insertion interval. In seconds.'), ...
            arg({'ignored','IgnoredEvents'},{},[],'Ignored event types. This is the list of event types that may occur between the opening and closing events. If any other event appears between a pair of successive opening and closing events, this range will not be considered for event insertion (it is considered "broken"). If set to ''ignoreall'', any event type may appear in between.', 'type','cellstr','shape','row')} ...
        },'Insertion interval definition. Events can be inserted either in fixed, absolute time window of the data set (absoluterange), or in time windows relative to reference events of a certain type (relativerange), or in time window spanned by two subsequent events of certain types (called the opening event and the closing event), optionally with ignored events in between (spannedrange).','cat','Time Ranges','mapper',@parse_segment), ...
    arg({'limits','Limits'},[-Inf Inf],[],'Time limits for event placement. Events that fall outside these bounds will be skipped. Therefore, intervals that intersect these boundaries may have fewer events than others.','cat','Time Ranges'),...
    arg({'event','InsertedType','Event'},'newevent',[],'Type of inserted events. The event type for the newly inserted events.','cat','Placement'), ...
    arg({'count','Count'},1,[],'Number of inserted events. This is either per interval or per second, depending on the Counting argument.','cat','Placement'), ...
    arg({'counting','Counting'},'perinterval',{'perinterval','persecond'},'Counting measure. Events can be inserted in a certain number per interval or per second','cat','Placement'), ...
    arg({'placement','Placement'},'equidistant',{'equidistant','random'},'Event placement scheme. Events can be inserted at equal (regular) distance from each other, or at random positions.','cat','Placement'), ...
    arg({'repeatable','Repeatable'},1,[],'Repeatable versus random placement. If 0, placement is random, if different from 0, the number is taken as the random seed, giving a unique repeatable run per number.','cat','Placement'), ...
    arg({'minlen','MinLength'},0,[],'Minimum segment length. Ignore segments that are shorter than this, in seconds.'), ...
    arg({'maxlen','MaxLength'},Inf,[],'Maximum segment length. Ignore segments that are longer than this, in seconds.'));

signal = opts.signal;
if isfield(signal,'epoch') && ~isempty(signal.epoch)
    error('the data set appears to contain epochs: only continuous data set are supported for now.'); end

% refine options
opts.limits = sort(max(min(opts.limits,signal.xmax),signal.xmin));
opts.limits = opts.limits*signal.srate;
opts.segmentspec.lo = opts.segmentspec.lo*signal.srate;
opts.segmentspec.hi = opts.segmentspec.hi*signal.srate;

% init randomization
if strcmp(opts.placement,'random') && opts.repeatable
    if hlp_matlab_version < 707
        % save & override RNG state
        randstate = rand('state'); %#ok<RAND>
        rand('state',5182+opts.repeatable); %#ok<RAND>
    else
        % create a legacy-compatible RandStream
        tracking.temp.randstream_inject_events = RandStream('swb2712','Seed',5182+opts.repeatable);
    end
end

newsignal = signal;
switch opts.segmentspec.arg_selection 
    case 'absoluterange'
        if opts.segmentspec.lo == -Inf
            opts.segmentspec.lo = signal.xmin*signal.srate; end
        if opts.segmentspec.hi == Inf
            opts.segmentspec.hi = signal.xmax*signal.srate; end
        if opts.segmentspec.hi < 0
            opts.segmentspec.hi = signal.xmax*signal.srate - opts.segmentspec.hi; end
        
        % inject using absolute latencies
        if isempty(opts.event)
            error('an event type must be specified'); end
        if ischar(opts.event)
            opts.event = make_default_event(signal.event,opts.event); end
        newsignal = perform_injection(newsignal,[opts.segmentspec.lo opts.segmentspec.hi],opts);
    case 'relativerange'
        % inject relative to one single marker
        if isempty(opts.event)
            opts.event = [opts.segmentspec.event '_i']; end
        if ischar(opts.event)
            opts.event = make_default_event(signal.event,opts.event); end
        for e=1:length(signal.event)
            if strcmp(signal.event(e).type,opts.segmentspec.event)
                newsignal = perform_injection(newsignal,signal.event(e).latency+[opts.segmentspec.lo opts.segmentspec.hi],opts); end
        end
	case 'spannedrange'
        % inject in between two successive markersputting that
        if isempty(opts.event)
            opts.event = [opts.segmentspec.openevent '_' opts.segmentspec.closeevent]; end
        if ischar(opts.event)
            opts.event = make_default_event(signal.event,opts.event); end
        startlat = [];
        for e=1:length(signal.event)
            % ignorance turned off but intermediate marker found?
            if ~iscell(opts.segmentspec.ignored) && ~strcmp(signal.event(e).type,opts.segmentspec.closeevent)
                startlat = []; end
            % selective ignorance turned on but non-ignored intermediate marker found?
            if ~isequal(opts.segmentspec.ignored,{'ignoreall'}) && ~any(strcmp(signal.event(e).type,[opts.segmentspec.ignored {opts.segmentspec.closeevent}]))
                startlat = []; end
            % found the begin-marker: open a new segmentspec
            if strcmp(signal.event(e).type,opts.segmentspec.openevent)
                startlat = signal.event(e).latency; end
            % reached the end-marker & have an open segmentspec: inject.
            if ~isempty(startlat) && strcmp(signal.event(e).type,opts.segmentspec.closeevent)
                newsignal = perform_injection(newsignal,[startlat+opts.segmentspec.lo,signal.event(e).latency+opts.segmentspec.hi],opts); end
        end
end
% sort the events by latency...
newsignal.event = newsignal.event(hlp_getresult(2,@sort,[newsignal.event.latency]));

% update .urevent field if trivial
if isempty(newsignal.urevent) || isequal([newsignal.event.urevent],1:length(newsignal.event))
    newsignal.urevent = newsignal.event;
    [newsignal.event.urevent] = arraydeal(1:length(newsignal.event));
end

% conclude randomization
if strcmp(opts.placement,'random') && opts.repeatable && hlp_matlab_version < 707
    % restore saved RNG state
    rand('state',randstate); %#ok<RAND>
end

exp_endfun;



function [signal,coverage] = perform_injection(signal,ival,opts)
global tracking;
coverage = 0;
% sanity check
if length(ival) == 2 && ival(1) <= ival(2)
    coverage = ival(2)-ival(1);
    % check the segment length
    seg_length = coverage/signal.srate;
    if seg_length < opts.minlen || seg_length > opts.maxlen
        return; end
    % hande the spacing method
    if strcmp(opts.counting,'persecond')
        opts.count = max(1,round(opts.count*(ival(2)-ival(1))/signal.srate)); end
    % handle the placement method
    if strcmp(opts.placement,'equidistant')
        if coverage ~= 0
            stepsize = (ival(2)-ival(1))/(opts.count-1);
            ival = round(ival(1):stepsize:ival(2));
        else
            ival = ival(1)*ones(1,opts.count);
        end
    elseif strcmp(opts.placement,'random')        
        if ival(1) < ival(2)            
            if hlp_matlab_version < 707
                positions = rand(1,opts.count);
            else
                positions = rand(tracking.temp.randstream_inject_events,1,opts.count);
            end
            ival = round(positions*(ival(2)-ival(1))+ival(1));
        elseif ival(1) == ival(2)
            ival = ival(1)*ones(1,opts.count);
        else
            ival = [];
        end
    else
        error('unsupported placement scheme specified');
    end
    if ~isempty(ival)
        % compute the individual latencies
        lats = ival(ival>=opts.limits(1) & ival<=opts.limits(2));
        % sanitize latencies
        lats = min(max(lats,1),signal.pnts);
        if ~isempty(signal.urevent)
            signal.urevent = []; end
        if isempty(signal.event)
            signal.event = setfield(setfield(opts.event,'latency',1),'type','dummy'); end; %#ok<SFLD>
        range = length(signal.event) + (1:length(lats));
        [signal.event(range)] = deal(opts.event);
        [signal.event(range).latency] = arraydeal(lats);
    end
end

% create a default event from an event array
function evt = make_default_event(evts,type)
if isempty(evts)
    evt = struct('type',{type},'latency',{[]},'duration',{1},'urevent',{[]});
else
    evt = evts(1);
    for fn=fieldnames(evt)'
        evt.(fn{1}) = []; end
    evt.type = type;
    if isfield(evt,'duration')
        evt.duration = 1; end
end


% parse a segmentspec specification into a cell array {tag,name,value,name,value,...}
function [selection,spec] = parse_segment(spec)
if isempty(spec)
    selection = 'absoluterange';
elseif ischar(spec) && any(strcmp(spec,{'absoluterange','relativerange','spannedrange'}))
	selection = spec;
    spec = {};
elseif iscell(spec) && any(strcmp(spec{1},{'absoluterange','relativerange','spannedrange'})) 
    selection = spec{1};
    spec = {};
elseif isfield(spec{1},'arg_selection')
    selection = spec{1}.arg_selection;
elseif any(strcmp(spec{1},{'absoluterange','relativerange','spannedrange'}))
    [selection,spec] = deal(spec{1},spec(2:end));
else
    % we have a custom segmentspec specification (as indicated in the function's help text)
    % parse it.
    mrks = {};
    lats = [];
    ignored = [];
    for i=1:length(spec)
        if ischar(spec{i})
            mrks{end+1} = spec{i};
        elseif iscell(spec{i})
            ignored = [ignored spec{i}];
        else
            lats(end+1) = spec{i};
        end
    end
    % and configure parameters
    if isempty(mrks)
        selection = 'absoluterange'; spec = {'lo' min(lats) 'hi' max(lats)};
    elseif length(mrks) == 1
        selection = 'relativerange'; spec = {'event' mrks{1} 'lo' min(lats) 'hi' max(lats)};
    elseif length(mrks) == 2
        if isequal(ignored,{})
            ignored = {'ignoreall'};
        elseif isequal(ignored,[])
            ignored = {};
        end
        selection = 'spannedrange'; spec = {'openevent' mrks{1} 'closeevent' mrks{2} 'lo' lats(1) 'hi' lats(2) 'ignored' ignored};
    end
end