function signal = set_new(varargin)
% Create a new EEGLAB data set from individual fields.
% Dataset = set_new(Arguments...)
%
% In:
%   Fields  : Pairs of field names and field values to add to the data set. fields not specified are
%             taken from eeg_emptyset, later fields override earlier fields; giving a struct in
%             place of a name-value pair is equivalent to writing out all the struct fieldnames and
%             respective values. fields that can be derived from others are derived.
%
%   `         optional special semantics:
%             * 'chanlocs' can be specified as cell-string array, and is generally completed using a 
%                default lookup
%             * 'data' can be specified as a cell array of data arrays, then concatenated across 
%                epochs, and with .epoch.target derived as the index of the cell which contained the 
%                epoch in question.
%             * 'tracking.online_expression' can be specified to override the online processing 
%                description
%
% Out:
%   Dataset : newly created EEGLAB set
% 
% Example:
%   % create a new continuous data set (with channels A, B, and C, and 1000 Hz sampling rate)
%   myset = set_new('data',randn(3,100000), 'srate',1000,'chanlocs',struct('labels',{'A','B','C'}));
%
%   % as before, but now also put in some events at certain latencies (note: latencies are in samples)
%   events = struct('type',{'X','Y','X','X','Y'},'latency',{1000,2300,5000,15000,17000});
%   myset = set_new('data',randn(3,100000), 'srate',1000, 'chanlocs',struct('labels',{'A','B','C'}), 'event',events);
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-28
dp;

% set_new_version<1.01> -- for the cache

if ~exp_beginfun('import') return; end

declare_properties('independent_channels',false,'independent_trials',false);

% construct the data set from arguments and defaults
signal = arg_define('allow-unlisted-names',varargin, ...
    ... % data description
    arg({'setname','DatasetName'},'new set',[],'Name of the new dataset.'), ...
    arg({'filename','FileName'},'',[],'Name of the file on disk. This is the file name (excluding the path) under which the set has previously been saved/loaded, if any.'), ...
    arg({'filepath','FilePath'},'',[],'Corresponding file path on disk. This is the path (excluding file name) under which the set has previously been saved/loaded, if any.'), ...
    arg({'subject','SubjectName'},'',[],'Subject name or identifier. This is the name or identifier of the person from whom the data was recorded.'), ...
    arg({'group','GroupName'},'',[],'Group name or identifier. This is the name or identifier of the group to which the subject belongs.'), ...
    arg({'condition','ConditionName'},'',[],'Condition name or identifier. This is the name or identifier of the experimental condition under which the data was recorded.','typecheck',false), ...
    arg({'session','SessionIdentifier'},[],[],'Identifier or number of the experiment session.','typecheck',false), ...
    arg({'comments','DatasetComments'},'created by set_new()',[],'Comments about the data set.'), ...
    ... % redundant data size fields (auto-deduced)
    arg_nogui({'nbchan','ChannelCount'},0,uint32([0 1000000]),'Number of channels in the data. Auto-deduced if omitted.'), ...
    arg_nogui({'trials','TrialCount'},0,uint32([0 1000000]),'Number of trials (epochs) in the data. Auto-deduced if omitted.'), ...
    arg_nogui({'pnts','SampleCount'},0,uint32([0 1000000000]),'Number of samples (time points) in the data. Auto-deduced if omitted.'), ...
    ... % time axis description
    arg({'srate','SamplingRate'},0,[0 0.1 100000 10000000],'Sampling rate in Hz.'), ...
    arg({'xmin','FirstTimepoint'},0,[],'First time-point in seconds. This is the time stamp of the first sample.'), ...
    arg_nogui({'xmax','LastTimepoint'},0,[],'Last time-point in seconds. This is the time stamp of the last sample. If given, should be consistent with xmin, srate and the size of the data.'), ...
    arg_nogui({'times','Timepoints'},[],[],'Time-points of the data. A vector of per-sample time points. Can be omitted.'), ...
    ... % signal payload
    arg({'data','DataArray'},[],[],'The matrix or tensor of data. Note: in previous versions a cell array was allowed here, too; if you need this functionality, set the typecheck flag of this argument to false.','shape','tensor'), ...
    ... % component-related fields
    arg_nogui({'icaact','ComponentActivationArray'},[],[],'A matrix/tensor of component activations. Typically requires that the other component-related fields are also specified.','shape','tensor'), ...
    arg_nogui({'icawinv','ComponentInverseWeights'},[],[],'Matrix of component spatial filter inverses. These are the forward projections of the filters. Typically equivalent to inv(.icaweights*.icasphere). In case of an ICA, this is the "sphering matrix".'), ...
    arg_nogui({'icasphere','ComponentSpheringMatrix'},[],[],'Matrix of component sphering transform. This is a linear transform (left-multiplied) used to sphere (de-correlate) the data prior to application of ica weights.'), ...
    arg_nogui({'icaweights','ComponentWeightMatrix'},[],[],'Matrix of component spatial filters. This is a linear transform (left-multiplied) used to derive component activations from the data after the sphering matrix has been applied. In case of an ICA, this is the "unmixing matrix".'), ...
    arg_nogui({'icachansind','ComponentChannelIndices'},[],uint32([1 1000000]),'Channel subset to derive components from. These are the indices of the data channels that shall be retained prior to applying the sphering matrix.','shape','row'), ...
    ... % channel descriptions
    arg({'chanlocs','ChannelLocations'},[],[],'Struct array of channel locations. May also be a cell array of channel labels.','type','expression'), ...
    arg_nogui({'urchanlocs','OriginalChannelLocations'},[],[],'Struct array of original channel locations. Describes the state of channel locations prior to processing.','type','expression'), ...
    arg_sub({'chaninfo','ChannelInfo'},{},{ ...
        arg_deprecated({'plotrad','PlottingRadius'},[],[],'Preferred plotting radius for topoplot().'), ...
        arg_deprecated({'shrink','ShrinkFactor'},[],[],'Preferred shrink factor for topoplot().'), ...
        arg({'nosedir','NoseDirection'},'+X',{'+X','+Y','+Z','-X','-Y','-Z'},'Direction in which the nose points. For correct orientation of the channel coordinate system.'), ...
        arg_deprecated({'nodatchans','NoDataChans'},[],[],'Legacy option for EEGLAB.'), ...
        arg_deprecated({'icachansind','ComponentChannelIndices'},[],uint32([1 1000000]),'Secondary location for component channel subset. Usage not recommended.'), ...
        arg({'labelscheme','LabelingScheme'},'',[],'Channel labeling scheme. For instance, ''10-20'' to indicate the 10-20 labeling system.'), ...
    },'Channel coordinate system information.','fmt','allow-unlisted-names'), ...
    arg_nogui({'ref','Reference'},'',[],'Reference channel. Can be ''common'', ''nasion'', a channel label, etc. Not used by BCILAB.'), ...
    ... % event descriptions
    arg({'event','Events'},[],[],'Event structure. This is a struct array with mandatory fields .type (string) and .latency (samples) per event.','typecheck',false,'shape','row'), ...
    arg_nogui({'urevent','OriginalEvents'},[],[],'Original event structure. Describes the state of events prior to processing.','typecheck',false,'shape','row'), ...
    arg_nogui({'eventdescription','EventDescription'},{{}},[],'Event description. Legacy EEGLAB option.','type','expression'), ...
    ... % epoch descriptions
    arg({'epoch','Epochs'},[],[],'Epoch structure. This is a struct array with information per epoch.','typecheck',false,'shape','row'), ...
    arg_nogui({'epochdescription','EpochDescription'},{{}},[],'Epoch description. Legacy EEGLAB option.','type','expression'), ...
    ... % miscellaneous EEGLAB fields
    arg_nogui({'reject','Rejections'},[],[],'Rejection information. Contains information about what data was rejected. Legacy.','type','expression'), ...
    arg_nogui({'stats','Statistics'},[],[],'Statistics information. Legacy.','type','expression'), ...
    arg_nogui({'specdata','SpecData'},[],[],'Spectral activation data. Legacy.'), ...
    arg_nogui({'specicaact','SpectralComponentActivations'},[],[],'Spectral component activation data. Legacy.'), ...
    arg_nogui({'splinefile','SplineFile'},'',[],'Cached spline file for plotting. Used by EEGLABs pop_headplot().'), ...
    arg_nogui({'icasplinefile','ComponentSplineFile'},'',[],'Cached spline file for plotting. Used by EEGLABs pop_headplot().'), ...
    ... % dipole fits
    arg({'dipfit','DipoleFits'},[],[],'Dipole fits per component.','typecheck',false,'shape','row'), ...
    ... % provenance tracking
    arg_nogui({'history','History'},'',[],'EEGLAB processing history script. Contains script commands that reproduce the data set.'), ...
    arg_nogui({'saved','WasSaved'},'no',{'no','yes','justloaded'},'Whether the data has been saved.'), ...
    ... % extension fields
    arg({'etc','Miscellaneous'},[],[],'Miscellaneous fields.','type','expression'));

% rewrite cell array data
if iscell(signal.data)
    if isempty(signal.data)
        signal.data = [];
    else
        sizes = cellfun('size',signal.data,1);
        if any(sizes(1) ~= sizes)
            error('For cell-array data, each cell must have the same number of channels.'); end
        sizes = cellfun('size',signal.data,2);
        if any(sizes(1) ~= sizes)
            error('For cell-array data, each cell must have the same number of time points.'); end
        data = signal.data;
        signal.data = [];
        for c=1:length(data)
            if isnumeric(data{c}) && ~isempty(data{c})
                signal.data = cat(3, signal.data, data{c});
                for i=size(signal.data,3)-size(data{c},3)+1:size(signal.data,3)
                    signal.epoch(i).target = c; end
            end
        end
    end
end

% process .chanlocs
if ~isfield(signal,'chanlocs') || isempty(signal.chanlocs)
    % create chanlocs from scratch, according to the data size
    if ~isempty(signal.data)
        signal.chanlocs = struct('labels',cellfun(@num2str,num2cell(1:size(signal.data,1),1),'UniformOutput',false),'type',repmat({'unknown'},1,size(signal.data,1))); 
    else
        signal.chanlocs = struct('labels',{},'type',{});
    end
else
    % bring chanlocs into an appropriate format
    try 
        signal.chanlocs = hlp_microcache('set_new1',@set_infer_chanlocs,signal.chanlocs); 
    catch e
        error('Could not look up channel locations according to the given chanlocs argument with error: %s (chanlocs were: %s)',e.message,hlp_tostring(signal.chanlocs,10000));
    end
end

% derive .xmax, .nbchan, .pnts, .trials
[signal.nbchan,signal.pnts, signal.trials, extra_dims] = size(signal.data); %#ok<NASGU>
signal.xmax = signal.xmin + (signal.pnts-1)/signal.srate;

% if epoched and there are events, derive the .epoch field
if signal.trials > 1 && ~isempty(signal.event) && isempty(signal.epoch)
    try
        signal = eeg_checkset(signal,'eventconsistency'); 
    catch e
        disp_once('set_new(): could not derive .epoch field due to error: %s',hlp_handleerror(e));
    end
end

% add .epoch.latency if possible
if ~isfield(signal.epoch,'latency')
    for i=1:length(signal.epoch)
        try
            tle = [signal.epoch(i).eventlatency{:}]==0;
            if any(tle)
                signal.epoch(i).latency = b.event(b.epoch(i).event(tle)).latency; end
        catch
        end
    end
end

% create .urevent field if applicable
if isempty(signal.urevent) && ~isempty(signal.event)
    signal.urevent = signal.event;
    [signal.event.urevent] = arraydeal(1:length(signal.event));
end

% do minimal consistency checks
if ~isempty(signal.chanlocs) && ~isempty(signal.data) && (length(signal.chanlocs) ~= signal.nbchan)
    error('The number of supplied channel locations (%i) does not match the number of channels (%i) in the data.',length(signal.chanlocs),signal.nbchan);  end
if isfield(signal,'epoch') && ~isempty(signal.epoch) && length(signal.epoch) ~= size(signal.data,3)
    error('The number of data epochs (%i) does not match the number of entries in the epoch field (%i).',size(signal.data,3),length(signal.epoch)); end

% ensure that .tracking.timeseries_fields is present
if ~isfield(signal,'tracking') || ~isfield(signal.tracking,'timeseries_fields')
    signal.tracking.timeseries_fields = {}; end

exp_endfun;