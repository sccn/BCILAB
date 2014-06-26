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
%                     can use fractional placement, which is based on the sampling rate  (most
%                     precise, but requires that the sampling rate is well within 1% of the true
%                     value) or placement next to the nearest sample (also works for streams that
%                     have irregular sampling rate). (default: 'fractional')
%
%   ClockAlignment : Clock alignment algorithm. The algorithm used to smooth clock alignment
%                    measurements; if the clock drift between data and markers is assumed to be low
%                    (e.g., come from same machine), then median is the safest choice. If drift is
%                    substantial (e.g., on a networked installation) then one may use linear to
%                    correct for that, but if the network is under heavy load it is safer to use
%                    the trimmed or robust estimators. The trimmed estimator survives occasional
%                    extreme load spikes, and the robust estimator further improves that tolerance
%                    by 2x. The only issue with the robust estimator is that whenever it updates
%                    (every 5s) it delays the BCI output by an extra 15ms of processing time.
%                    (default: 'trimmed')
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
        arg({'update_freq','UpdateFrequency'},20,[0 Inf],'Update frequency. The rate at which new chunks of data is polled from the device, in Hz.'), ...
        arg({'buffer_len','BufferLength'},10,[],'Internal buffering length. This is the maximum amount of backlog that you can get.'), ...
        arg({'channel_override','ChannelOverride'}, [], [], 'Override channel labels. This allows to replace the channel labels that are provided by the stream.','type','cellstr','shape','row'), ...
        arg({'srate_override','SamplingRateOverride'}, 0, [], 'Override sampling rate. This allows to replace the sampling rate that is provided by the stream.'), ...
        arg({'marker_placement','MarkerPlacement'}, 'nearest', {'nearest','interpolated'}, 'Marker placement rule. Controls how the latency of markers is determined -- can use interpolated placement, which is based on the sampling rate (but requires that the sampling rate is well within 1% of the true value) or placement next to the nearest sample (works for streams that have irregular sampling rate but fails for stream which incorrectly report their sampling rate).','guru',true), ...
        arg({'clock_alignment','ClockAlignment'}, 'median', {'trimmed','median','robust','linear','raw'}, 'Clock alignment algorithm. The algorithm used to smooth clock alignment measurements; if the clock drift between data and markers is assumed to be low (e.g., come from same machine), then median is the safest choice. If drift is substantial (e.g., on a networked installation) then one may use linear to correct for that, but if the network is under heavy load it is better to use the trimmed or robust estimator. The trimmed estimator survives occasional extreme load spikes, and the robust estimator further improves that tolerance by 2x. The only issue with the robust estimator is that whenever it updates (every 5s) it delays the BCI output by an extra 15ms of processing time. Most estimators other than median can introduce a few-ms jitter between data and markers.','guru',true), ...
        arg({'jitter_correction','JitterCorrection'}, true, [], 'Correct for jittered time stamps. This applies only to the data stream.','guru',true), ...
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
    
    % create a new inlet
    disp('Opening an inlet...');
    inlet = lsl_inlet(result{1});
    % get the stream info...
    info = inlet.info();

    % check srate
    if info.nominal_srate() <= 0
        warning('BCILAB:run_readlsl:not_supported','The given stream has a variable sampling rate. This plugin currently only supports data with regular sampling rate.'); end

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
        if nominal_srate && strcmp(opts.marker_placement,'nearest')
            fprintf('IMPORTANT: If you know that your stream''s sampling rate is specified incorrectly and override it, it is recommended to also set MarkerPlacement to ''fractional'', which bases it on your given sampling rate (nearest will most likely give incorrect results).'); end
        nominal_srate = opts.srate_override; 
    end

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

    % initialize shared variables
    samples = {};
    timestamps = [];
    
    % start background acquisition on the online stream (set up read_data as callback function)
    onl_read_background(opts.new_stream,@()read_data(inlet,marker_inlet,opts.always_double),opts.update_freq);

    disp('Now reading...');
    
    
    function result = read_data(data_inlet,marker_inlet,always_double)
        % get a new chunk of data
        [chunk,stamps] = data_inlet.pull_chunk();
        stamps = stamps + data_inlet.time_correction([],opts.clock_alignment);
        if always_double
            chunk = double(chunk); end
        
        if ~isempty(marker_inlet) && ~isempty(stamps)
            % receive any available markers
            corr = marker_inlet.time_correction([],opts.clock_alignment);
            while true
                [sample,ts] = marker_inlet.pull_sample(0.0);
                if ts                
                    samples(end+1) = sample; %#ok<AGROW>
                    timestamps(end+1) = ts+corr; %#ok<AGROW>
                else
                    break;
                end
            end

            % submit all markers that overlap the chunk being submitted (some may be ahead of the data stream received thus far)
            matching = timestamps < stamps(end);
            if any(matching)
                % calculate marker latencies, in samples relative to the beginning of the chunk being submitted
                if (strcmp(opts.marker_placement,'nearest') || ~nominal_srate)
                    latencies = argmin(abs(bsxfun(@minus,timestamps(matching),stamps(:))));
                else
                    latencies = 1 + (timestamps(matching)-stamps(1))*nominal_srate;
                end
                latencies = min(size(chunk,2),max(1,latencies));
                markers = struct('type',samples(matching),'latency',num2cell(latencies));
                % and remove the markers from the buffer
                samples(matching) = [];
                timestamps(matching) = [];
            else
                markers = [];
            end
            % construct output
            result = {chunk,markers};
        else
            result = chunk;
        end
    end

end
