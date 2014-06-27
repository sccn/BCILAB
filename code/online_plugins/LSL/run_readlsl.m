function run_readlsl(varargin)
% Receive real-time data from a source on the lab streaming layer.
% run_readlsl(MatlabStream,SelectionProperty,SelectionValue,UpdateFrequency)
%
% This plugin connects to and receives data from a device on the lab streaming layer. The data source
% is specified by means of a query (some allowed properties are type, name, channel_count, and srate).
%
% In:
%   StreamName : name of the stream; a variable with this name will be created in the MATLAB workspace 
%                to hold the stream's data. If such a variable already exists it will be overridden.
%
%   DataStreamQuery : Data stream query. Allows to select an LSL data stream to read from (e.g., by setting 
%                     it to 'type=''EEG''' or 'name=''BioSemi'''). (default: 'type=''EEG''')
%
%   MarkerStreamQuery : Marker stream query. Allows to select an LSL marker stream to read from (e.g., by setting 
%                       it to 'type=''Markers'''). Leave it empty to ignore markers. (default: 'type=''Markers''')
%
%   ConvertToDouble : Always convert the signal to double precision. (default: true)
%   
%   UpdateFrequency : The rate at which new chunks of data is polled from the device, in Hz. 
%                     (default: 20)
%
%   BufferLength : Internal buffering length. This is the maximum amount of backlog that you can
%                  get, in seconds. (default: 30)
%
%   ChannelOverride : Override channel labels. This allows to replace the channel labels that 
%                     are provided by the stream. (default: {})
%
%   SamplingRateOverride : Override sampling rate. This allows to replace the sampling rate that is
%                          provided by the stream. (default: 0)
%
%   MarkerPlacement : Marker placement rule. Controls how the latency of markers is determined --
%                     can use fractional placement, which is based on the sampling rate (potentially
%                     the most precise, but requires that the sampling rate is well within 1% of the
%                     true value) or placement next to the nearest sample (also works for streams
%                     that have irregular sampling rate). (default: 'nearest')
%
%   ClockAlignment : Clock alignment algorithm. The algorithm used to smooth clock alignment
%                    measurements; if the clock drift between data and markers is assumed to be low
%                    (e.g., come from same machine), then median is the safest choice. If drift is
%                    substantial (e.g., on a networked installation) then one may use linear to
%                    correct for that, but if the network is under heavy load it is safer to use
%                    the trimmed or robust estimators. The trimmed estimator survives occasional
%                    extreme load spikes, and the robust estimator further improves that tolerance
%                    by 2x. The only issue with the robust estimator is that whenever it updates
%                    (every 5s) it delays the BCI output by an extra 15ms of processing time. Zero
%                    disables the time-correction. (default: 'median')
%
%   JitterCorrection : Correct for jittered time stamps. This corrects jitter in the time stamps of
%                      the data stream chunks assuming that the underlying sampling rate is regular
%                      but may drift slowly. (default: true)
%
%   ForgetHalftime : Forget factor as information half-life. In estimating the effective sampling rate 
%                    a sample which is this many seconds old will be weighted 1/2 as much as the
%                    current sample in an exponentially decaying window. (default: 30)
%
% Notes:
%   The general format of the queries is that of an XPath 1.0 predicate on the meta-data of a given stream.
%
%   Some older versions of MATLAB (solved since 2009a) cannot handle a sender disconnecting while
%   run_readlsl is running; on those versions you need to clear (i.e., stop) the stream before you
%   disconnect at the sender side, otherwise your MATLAB will hang.
%
% Examples:
%   % receive data from a device that contains gaze data
%   run_readlsl('mystream','type','Gaze');
%
%   % read from an EEG stream (default) but use custom channel labels
%   run_readlsl('ChannelOverride', {'C3','C4','Cz','O1','O2'};
%
%   % read from a stream that has the name 'PhaseSpace'
%   run_readlsl('SelectionProperty','name','SelectionValue','PhaseSpace')
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2012-03-21

    persistent lib;

    % declare the name of this component (shown in the menu)
    declare_properties('name','Lab streaming layer');

    % read options
    opts = arg_define(varargin, ...
        arg({'new_stream','MatlabStreamName','MatlabStream'}, 'laststream',[],'MATLAB Stream Name. A variable with this name will be created in the MATLAB workspace to hold the stream''s data. If such a variable already exists it will be overridden.','type','char'), ...
        arg({'data_query','DataStreamQuery','DataQuery'}, 'type=''EEG''',[],'Data stream query. Allows to select a data stream to read from (e.g., by setting it to type=''EEG'' or name=''BioSemi'').'), ...
        arg({'marker_query','MarkerStreamQuery','MarkerQuery'},'type=''Markers''',[],'Marker stream query. Allows to select a marker stream to read from (e.g., by setting it to type=''Markers''). Leave it empty to ignore markers.','shape','row','type','char'), ...
        arg({'always_double','ConvertToDouble'},true,[],'Convert to double. Always convert the signal to double precision.'), ...
        arg({'update_freq','UpdateFrequency'},20,[0 0.0001 200 Inf],'Update frequency. The rate at which new chunks of data is polled from the device, in Hz.'), ...
        arg({'buffer_len','BufferLength'},10,[0 0.1 300 Inf],'Internal buffering length. This is the maximum amount of backlog that you can get.'), ...
        arg({'channel_override','ChannelOverride'}, [], [], 'Override channel labels. This allows to replace the channel labels that are provided by the stream.','type','cellstr','shape','row'), ...
        arg({'srate_override','SamplingRateOverride'}, 0, [], 'Override sampling rate. This allows to replace the sampling rate that is provided by the stream.'), ...
        arg({'marker_placement','MarkerPlacement'}, 'nearest', {'nearest','interpolated'}, 'Marker placement rule. Controls how the latency of markers is determined -- can use interpolated placement, which is based on the sampling rate (but requires that the sampling rate is well within 1% of the true value) or placement next to the nearest sample (works for streams that have irregular sampling rate).','guru',true), ...
        arg({'clock_alignment','ClockAlignment'}, 'median', {'trimmed','median','robust','linear','raw','zero'}, 'Clock alignment algorithm. The algorithm used to smooth clock alignment measurements; if the clock drift between data and markers is assumed to be low (e.g., come from same machine), then median is the safest choice. If drift is substantial (e.g., on a networked installation) then one may use linear to correct for that, but if the network is under heavy load it is better to use the trimmed or robust estimator. The trimmed estimator survives occasional extreme load spikes, and the robust estimator further improves that tolerance by 2x. The only issue with the robust estimator is that whenever it updates (every 5s) it delays the BCI output by an extra 15ms of processing time. Most estimators other than median can introduce a few-ms jitter between data and markers. Zero disables the time-correction.','guru',true), ...
        arg({'jitter_correction','JitterCorrection'}, true, [], 'Correct for jittered time stamps. This corrects jitter in the time stamps of the data stream chunks assuming that the underlying sampling rate is regular.','guru',true), ...
        arg({'forget_halftime','ForgetHalftime'}, 30, [1 20 60 Inf], 'Forget factor as information half-life. In estimating the effective sampling rate a sample which is this many seconds old will be weighted 1/2 as much as the current sample in an exponentially decaying window.','guru',true), ...
        arg({'clock_est_window','ClockEstimationWindow'}, 30, [1 10 120 Inf], 'Clock alignment estimation window. This is the time window over which clock alignment between the marker stream and data stream will be estimated (using the estimator selected by the ClockAlignment option). Note that JitterCorrection, when enabled, will implicitly also smooth any jitter in the clock alignment.','guru',true), ...
        arg_deprecated({'property','SelectionProperty'}, '',[],'Selection property. The selection criterion by which the desired device is identified on the net. This is a property that the desired device must have (e.g., name, type, desc/manufacturer, etc.'), ...
        arg_deprecated({'value','SelectionValue'}, '',[],'Selection value. This is the value that the desired device must have for the selected property (e.g., EEG if searching by type, or Biosemi if searching by manufacturer).'));

    % get a library handle (here with an explicit path because we want it to work if the toolbox is compiled, too)
    if isempty(lib)
        lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin')); end

    if ~isempty(opts.property) && ~isempty(opts.value)
        opts.data_query = [opts.property '=''' opts.value '''']; end

    % look for the desired device
    disp(['Looking for a device with ' opts.data_query ' ...']);
    result = {};
    while isempty(result)
        result = lsl_resolve_bypred(lib,opts.data_query); end
    
    % create a new inlet & query stream info
    disp('Opening an inlet...');
    inlet = lsl_inlet(result{1});
    info = inlet.info();

    % try to get the channel labels & sanity-check them
    if isempty(opts.channel_override)
        channels = {};    
        ch = info.desc().child('channels').child('channel');
        if ch.empty()
            ch = info.desc().child('channel'); end
        while ~ch.empty()
            name = ch.child_value_n('label');
            if isempty(name)
                name = ch.child_value_n('name'); end
            if name
                channels{end+1} = name; end %#ok<AGROW>
            ch = ch.next_sibling_n('channel');
        end
    else
        channels = opts.channel_override;
    end
    if length(channels) ~= info.channel_count()
        disp('The number of channels in the stream does not match the number of labeled channel records. Using numbered labels.');
        channels = cellfun(@(k)['Ch' num2str(k)],num2cell(1:info.channel_count(),1),'UniformOutput',false);
    end

    % allow the nominal_srate to be overridden
    nominal_srate = info.nominal_srate();
    if ~nominal_srate && ~opts.srate_override
        fprintf('The given data stream reports a sampling rate of 0 (= irregular sampling); if you know that the stream is actually regularly sampled we recommend that you override the sampling rate by passing in a nonzero SamplingRateOverride value.'); end
    if opts.srate_override
        nominal_srate = opts.srate_override; end

    % create online stream data structure in base workspace (using appropriate meta-data)
    onl_newstream(opts.new_stream,'srate',nominal_srate,'chanlocs',channels,'buffer_len',opts.buffer_len);

    if ~isempty(opts.marker_query)
        % also try to find the desired marker stream
        disp(['Looking for a marker stream with ' opts.marker_query ' ...']);
        result = {};
        while isempty(result)
            result = lsl_resolve_bypred(lib,opts.marker_query); end
        % create a new inlet
        disp('Opening a marker inlet...');
        marker_inlet = lsl_inlet(result{1});
    else
        marker_inlet = [];
    end
    
    % initialize marker buffer
    marker_data = {};       % marker samples gathered so far and to be committed with next overlapping chunk
    marker_stamps = [];     % time stamps associated with those marker samples
    
    % state variables for recursive least squares jitter correction
    P = 1e10*eye(2);        % precision matrix (inverse covariance matrix of predictors)
    w = [0 0]';             % linear regression coefficients [offset,slope]
    lam = 2^(-1/(nominal_srate*opts.forget_halftime)); % forget factor in RLS calculation
    n = 0;                  % number of samples observed so far    
    numeric_offset = [];    % time-stamp offset to keep numerics healthy; will be initialized with first measured time stamp
    
    % start background acquisition on the online stream (set up read_data as callback function)
    onl_read_background(opts.new_stream,@()read_data(inlet,marker_inlet,opts.always_double),opts.update_freq);

    disp('Now reading...');
    
    
    function result = read_data(data_inlet,marker_inlet,always_double)
        % get a new chunk of data
        [chunk,stamps] = data_inlet.pull_chunk();
        data_clock = data_inlet.time_correction([],opts.clock_alignment,opts.clock_est_window);
        stamps = stamps + data_clock;
        if opts.jitter_correction
            stamps = update_regression(stamps); end
        if always_double
            chunk = double(chunk); end
        
        if ~isempty(marker_inlet) && ~isempty(stamps)
            % receive any available markers
            marker_clock = marker_inlet.time_correction([],opts.clock_alignment,opts.clock_est_window);
            while true
                [sample,ts] = marker_inlet.pull_sample(0.0);
                if ts                
                    marker_data(end+1) = sample; %#ok<AGROW>
                    marker_stamps(end+1) = ts + marker_clock; %#ok<AGROW>
                else
                    break;
                end
            end

            % submit all markers that overlap the current chunk (some may be ahead of the chunks
            % received thus far and will therefore be submitted on a subsequent update)
            matching = marker_stamps < stamps(end);
            if any(matching)
                % calculate marker latencies, in samples relative to the beginning of the chunk being submitted
                if (strcmp(opts.marker_placement,'nearest') || ~nominal_srate)
                    latencies = argmin(abs(bsxfun(@minus,marker_stamps(matching),stamps(:))));
                else
                    latencies = 1 + (marker_stamps(matching)-stamps(1))*nominal_srate;
                end
                latencies = min(size(chunk,2),max(1,latencies));
                markers = struct('type',marker_data(matching),'latency',num2cell(latencies));
                % and remove the markers from the buffer
                marker_data(matching) = [];
                marker_stamps(matching) = [];
            else
                markers = [];
            end
            % construct output
            result = {chunk,markers};
        else
            result = chunk;
        end
    end

    % perform RLS block update of regression coefficients
    % this is a regression from sample index onto timestamp of the sample
    function y = update_regression(y)
        if ~isempty(y)
            % sanitize numerics (all done relative to the first observed time stamp)
            if isempty(numeric_offset)
                numeric_offset = y(1); end
            y = y - numeric_offset;        
            % define predictor matrix (bias, sample index)
            X = [ones(1,length(y)); n + (1:length(y))];
            n = n + length(y);            
            % apply updates...
            for t=1:length(y)
                u = X(:,t);
                d = y(t);
                pi = u'*P;
                gam = lam + pi*u;
                k = pi'/gam;
                al = d - w'*u;
                w = w + k*al;
                Pp = k*pi;
                P = (1/lam)*(P-Pp);
            end            
            % predict y
            y = w'*X + numeric_offset;
        end
    end

end
