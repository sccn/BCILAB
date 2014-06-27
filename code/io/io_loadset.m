function res = io_loadset(varargin)
% Load/import a data set from disk, across all formats supported by EEGLAB.
% Set = io_loadset(Filename, Options...)
%
% This function is used to import continuous/raw data for processing with BCILAB. Generally, the 
% data can be loaded in reduced form (e.g., to conserve memory or to exclude misc channels), and 
% a variety of file formats support special options. These and all other options can be passed as
% name-value pairs following the file name. 
%
% Some additional support exists for extraction of trigger channels (which has default settings, but
% which may have to be customized if the trigger channel format is unusual), as well as for tagging
% data with meta-information (for possible use in study- level processing).
%
% In:
%   Filename      :   name of the file; platform-independent path preferred.
%
%   Options...    :   --- parameters for data reductions at load time ---
%                     'channels'    : channel index subset to be loaded (memory-efficient for .vhdr,
%                                     .eeg, .bdf (if BioSig works), .ctf, .ds)
%                     'samplerange' : sample range to be read; [first_sample last_sample]
%                                     (memory-efficient for .vhdr, .raw (EGI), .cnt (Neuroscan, EEProbe))
%                     'timerange'   : time range to be read, in seconds; [begin_time end_time]
%                                     (memory-efficient for .cnt (Neuroscan, EEProbe), .bdf, .ctf, .ds)
%                     'types'       : type(s) of channels to retain; either a string or cell array
%                                     of strings (e.g., 'EEG')
%                     'subsampled'  : sub-sample the data to the given rate at load time; note -- 
%                                     this is not the same as having a resample filter in the
%                                     processing pipeline: the filter is more general (e.g., will
%                                     apply properly during online processing while this operation
%                                     changes the data before the bcilab pipeline sees it) (default: [])
%
%                     --- parameters for dataset annotation ---
%                     'setname'     : initial data set name (default: '')
%                     'subject'     : initial subject identifier (default: '')
%                     'group'       : initial group identifier (default: '')
%                     'condition'   : initial condition identifier (default: '')
%                     'session'     : initial session number (default: [])
%                     'comments'    : initial data set comments (default: '')
%
%                     --- parameters for post-processing ---
%                     'markerchannel' : if a marker channel is present, allows to customize settings 
%                                       for how events shall be derived from it. This is a cell
%                                       array of parameters to set_infer_markers(), ideally name-value 
%                                       pairs.
%                     'montage_disambiguation' : rule for handling ambiguous montages (i.e., the file's
%                                                channel labels match multiple montage files similarly well):
%                                                * 'coverage' : picks the one with maximum coverage
%                                                * 'first' : picks the first one in the list of candidates
%                                                             (out of those that have reasonable coverage)
%                                                * 'deduce' : uses correlation analysis to deduce the best fit
%                                                * filename : uses the montage with the given file name
%                                                (default: 'deduce')
%
%                     --- format-specific importing parameters ---
%                    .sna: 'gain', see pop_snapread
%                    .cnt: 'keystroke','memmapfile','scale','dataformat','blockread', 'triggerfile' see pop_loadcnt/readcnt
%                    .eeg: 'range_trials','range_typeeeg','range_response','format', see pop_loadeeg
%                    .ctf: 'trials', see pop_readctf/ctf_read
%                    .ds: 'trials', see pop_ctf_read
%                    .dat: 'mergeposition','concatruns','maxevents', see BCI2000import
%                    .bdf: see pop_biosig/pop_readbdf
%                    .raw: 'samplerate' : the sampling rate of the data
%                    .xdf: 'streamname','streamtype','effective_rate','exclude_markerstreams'; see eeg_load_xdf
%
% Out:
%   Set :   data set in EEGLAB format
%
% Notes:
%   Note that the output of this function is not the data set itself but rather a "proxy" (or handle)
%   to it, which will be resolved by the toolbox into the corresponding EEGLAB set at the moment when 
%   it is needed. This is for efficiency (because most of the time, the raw data is not further needed
%   after the required pre-processed versions of it have been computed once and cached). You can 
%   manually resolve it into the actual data by calling mydata = exp_eval(mydata).
%
% Examples:
%   % load a given data set (here: BrainProducts format)
%   raw = io_loadset('/data/projects/test.vhdr')
%
%   % as before, but using a platform-independent file path (recommended)
%   % this assumes that BCILAB's data path has been initialized to /data/projects in the startup
%   % configuration (e.g. via the line: data = '/data/projects';)
%   raw = io_loadset('data:/test.vhdr')
%
%   % like before, but this time restrict the channel range to the first 32
%   raw = io_loadset('data:/test.vhdr','channels',1:32)
%
%   % customize how the marker channel in the data should be identified or parsed
%   % (for example, some recordings may have more events than the default cutoff for what is
%   % considered a marker channel; here, this cutoff is raised to 50000)
%   raw  io_loadset('/data/sets/myrecording.bdf','markerchannel',{'MaxEvents',50000})
%
% See also:
%   pop_loadset, pop_loadbv, pop_readegi, pop_read_erpss, pop_fileio,
%   pop_snapread, pop_loadcnt, pop_loadeep, pop_loadeeg, pop_readbdf,
%   pop_biosig, pop_loadbva, pop_ctf_read, BCI2000import,
%   set_infer_markers, set_infer_chanlocs
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-03-28

% io_loadset_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% read options
allopts = arg_define([0 1],varargin,...
    arg_norep({'filename','Filename'},mandatory,[],'File name to load. Any format that is recognized by EEGLAB should be supported.','type','char','shape','row'), ...
    ... % data reductions
    arg({'channels','ChannelSubset','chans'},[],[],'Channel index subset. Restrict the loaded channels to a subset.'), ...
    arg({'samplerange','SampleRange','srange'},[],[],'Sample range. Restrict the loaded data to a sub-range, in samples.'), ...
    arg({'timerange','TimeRange','trange'},[],[],'Time range. Restrict the loaded data to a sub-range, in seconds.'), ...
    arg({'types','ChannelTypes','type'},[],[],'Channel type subset. Restrict the loaded channels to those with matching types; should be a cell-array of strings.','type','cellstr'), ...
    arg({'casttodouble','CastToDouble'},true,[],'Cast data to double-precision. Takes up more space but helps with filters that are not working properly when a single-precision time series is passed in.'), ...
    arg({'subsampled','Subsampled'},[],[],'Sub-sample the data at load time. This is not the same as having a resample filter in the processing pipeline: the filter is more general (e.g., will apply properly during online processing while this operation changes the data before the bcilab processing pipeline or the data curation script sees it).'), ...
    ... % marker channel handling
    arg_sub({'markerchannel','MarkerChannel'},{},@set_infer_markers,'Marker channel processing. Optional parameters to control the marker channel processing function.'), ...
    ... % chanlocs handling
    arg({'infer_chanlocs','InferChanlocs'},true,[],'Infer channel locations if necessary. This will look up locations from standard tables if missing and attempt to deduce the channel types.'), ...
    arg({'montage_disambiguation','MontageDisambiguation'},'deduce',[],'Rule for handling ambiguous montages. This situation occurs when the file''s channel labels match multiple montage files similarly well). Options are: ''coverage'': picks the one with maximum coverage, ''deduce'': uses correlation analysis to deduce the best fit, filename : uses the montage with the given file name'), ...
    ... % added annotations
    arg({'setname','DatasetName'},'',[],'Data set name. This is meta-information of potential use in study-level processing.'), ...
    arg({'subject','SubjectIdentifier'},'',[],'Subject identifier. This is meta-information of potential use in study-level processing.'), ...
    arg({'group','GroupIdentifier'},'',[],'Group identifier. This is meta-information of potential use in study-level processing.'), ...
    arg({'condition','ConditionIdentifier'},'',[],'Condition identifier. This is meta-information of potential use in study-level processing.'), ...
    arg({'session','SessionIdentifier'},'',[],'Session identifier. This is meta-information of potential use in study-level processing.'), ...
    arg({'comments','DataComments'},'',[],'Data comments. This is meta-information.'), ...
    ... % misc arguments (neither displayed nor assigned by default)
    arg_norep('gain',unassigned), arg_norep('keystroke',unassigned), arg_norep('memmapfile',unassigned), arg_norep('scale',unassigned), arg_norep('dataformat',unassigned), ...
    arg_norep('blockread',unassigned), arg_norep('triggerfile',unassigned), arg_norep('range_trials',unassigned), arg_norep('range_typeeeg',unassigned), arg_norep('range_response',unassigned), ...
    arg_norep('format',unassigned), arg_norep('trials',unassigned), arg_norep('mergeposition',unassigned), arg_norep('concatruns',unassigned), arg_norep('maxevents',unassigned), ...
    arg_norep('samplerate',unassigned), arg_norep('blockrange',unassigned),arg_norep('ref',unassigned),arg_norep('rmeventchan',unassigned), arg_norep('exclude_markerstreams',unassigned),...
    arg_norep('streamname',unassigned), arg_norep('streamtype',unassigned), arg_norep('effective_rate',unassigned), ...
    arg_deprecated('nofixups',false,[],'This parameter has been retired.'));

opts = rmfield(allopts,{'filename','types','casttodouble','setname','subject','group','condition','session','comments','markerchannel','subsampled','infer_chanlocs','montage_disambiguation','nofixups'});
filename = env_translatepath(allopts.filename);
[base,name,ext] = fileparts(filename);

% test if the file actually exists
if ~exist(filename,'file')
    error(['The file ' filename ' does not exist.']); end

% ... and whether it can be opened
f = fopen(filename,'r');
if f==-1
    error(['The file ' filename ' could not be opened; please check your file permissions.']);
else
    fclose(f);
end

% sanity check
if ~isempty(opts.samplerange) && ~isempty(opts.timerange)
    error('Please do not specify both a sample subset and time subset.'); end

% bitrate test function (for some formats)
wrong_bitrate = @(eeg) mad(mean(eeg.data(:,1:2:end),2)./mean(eeg.data(:,2:2:end),2),1) > 0.5;

disp(['io_loadset(): loading ' filename '...']);
try
    % if the extension is supported, dispatch to the specific loader / importer
    switch lower(ext)
        case '.set'
            % EEGLAB data set
            args = hlp_struct2varargin(opts,'suppress',{'channels','samplerange','timerange'});
            try
                res = pop_loadset('filepath',[base filesep], 'filename', [name ext], args{:});
            catch e
                % fall back to direct import attempt
                fprintf('pop_loadset failed with error: %s\n',e.message);
                fprintf('attempting direct import...\n');
                res = getfield(io_load(filename,'-mat'),'EEG');
                if ischar(res.data)
                    binfile = [base filesep res.data];
                    if ~exist(binfile,'file')
                        error('The associated raw-data file %s was not found.',binfile); end
                    f = fopen(binfile,'r','ieee-le');
                    if f==-1
                        error('The associated raw-data file %s does could not be opened. Please check your file permissions.',binfile); end
                    if strcmp(res.data(end-3:end),'.fdt')
                        res.data = fread(f,[res.nbchan Inf], 'float32');
                    elseif strcmp(res.data(end-3:end),'.dat')
                        res.data = fread(f,[res.trials*res.pnts EEG.nbchan], 'float32')';
                    end
                    fclose(f);
                end
            end
            if ~isfield(res,'tracking') || ~isfield(res.tracking,'online_expression')
                % it comes fresh from EEGLAB
                disp('The loaded EEGLAB set is lacking an online expression; assuming it contains unfiltered data.')
                disp('If it contains filtered data, however, BCI models derived from it will likely not be online-capable.');
            else
                % it came out of a curation script in BCILAB: check the fingerprint before we sign it off as freshly loaded
                if ~isfield(res.tracking,'fingerprint') || hlp_fingerprint(rmfield(res,'tracking')) ~= res.tracking.fingerprint
                    disp('The loaded data set has been edited manually; forgetting its old tracking information.');
                    res = rmfield(res,'tracking');
                end
            end
        case '.vhdr'
            % BrainProducts Vision Recorder file
            res = pop_loadbv([base filesep],[name ext],opts.samplerange,opts.channels);
            opts.channels = [];
            opts.samplerange = [];
            % remove non-standard event fields
            res.event = rmfield(res.event,intersect({'code','bvtime','channel'},fieldnames(res.event)));
        case '.raw'
            try
                % EGI continuous raw data
                res = pop_readegi(filename,opts.samplerange);
                opts.samplerange = [];
            catch
                try
                    disp('EGI .raw importer failed; falling back to ERPSS .raw importer.');
                    if ~isfield(opts,'samplerate')
                        disp('ERPSS .raw formats can only be loaded if you explicitly pass a ''samplerate'' parameter.'); end
                    % ... or ERPSS raw
                    res = pop_read_erpss(filename,opts.samplerate);
                catch
                    disp('ERPSS .cnt importer failed; falling back to FileIO Yokogawa importer.');
                    % .. or Yokogawa (and possibly other formats...)
                    args = hlp_struct2varargin(opts,'suppress',{'channels','timerange','samplerange'});
                    res = pop_fileio(filename,args{:});
                end
            end
        case '.sna'
            % Snapmaster SNA
            args = hlp_struct2varargin(opts,'restrict',{'gain'});
            res = pop_snapread(filename,args{:});
        case '.cnt'
            % Neuroscan CNT
            try
                args = hlp_struct2varargin(opts,'suppress',{'channels','samplerange','timerange'});
                if ~isempty(opts.samplerange)
                    args = [args {'sample1',opts.samplerange(1),'ldnsamples',opts.samplerange(2)-opts.samplerange(1)+1}]; end
                if ~isempty(opts.timerange)
                    args = [args {'t1',opts.timerange(1),'lddur',opts.timerange(2)-opts.timerange(1)}]; end
                res = pop_loadcnt(filename, args{:});
                % check if it should have been 32 bits...
                if wrong_bitrate(res) && ~isfield(opts,'dataformat')
                    disp('The data is likely 32 bits; re-loading.');
                    res = pop_loadcnt(filename, args{:}, 'dataformat','int32');
                end
                opts.samplerange = [];
                opts.timerange = [];
            catch
                disp('Neuroscan importer failed; falling back to ANT EEProbe importer.');
                % ... or ANT EEProbe CNT
                args = hlp_struct2varargin(opts,'restrict',{'triggerfile'});
                if ~isempty(opts.samplerange)
                    args = [args {'sample1',opts.samplerange(1),'sample2',opts.samplerange(2)}]; end
                if ~isempty(opts.timerange)
                    args = [args {'time1',opts.timerange(1),'time2',opts.timerange(2)}]; end
                res = pop_loadeep(filename,args{:});
                opts.samplerange = [];
                opts.timerange = [];
            end
        case '.avr'

            % ANT EEProbe average file
            try
                res = pop_loadeep_avg(filename);
            catch
                % .. or Megis / BESA; via FileIO
                args = hlp_struct2varargin(opts,'suppress',{'channels','timerange','samplerange'});
                res = pop_fileio(filename,args{:});
            end
        case '.eeg'
            % Neuroscan EEG
            try
                optseeg = hlp_varargin2struct(opts,'range_trials',[],'range_typeeeg',[],'range_response',[],'format','short');
                res = pop_loadeeg([name ext],[base filesep],opts.channels,optseeg.range_trials,optseeg.range_typeeg,optseeg.range_response,optseeg.format);
                if wrong_bitrate(res) && ~isfield(opts,'format')
                    disp('The data is likely 32 bits; re-loading.');
                    res = pop_loadeeg([name ext],[base filesep],opts.channels,optseeg.range_trials,optseeg.range_typeeg,optseeg.range_response,'int32');
                end
                opts.channels = [];
            catch
                % ... or ANT / BrainProducts file; via FileIO
                args = hlp_struct2varargin(opts,'suppress',{'channels','timerange','samplerange'});
                res = pop_fileio(filename,args{:});
            end
        case {'.bdf','.edf'}
            % BioSEMI BDF/EDF
            try
                args = hlp_struct2varargin(opts,'suppress',{'channels','samplerange','timerange'}, 'rewrite',{'timerange','range'});
                res = pop_readbdf(filename,args{:});
                opts.timerange = [];
            catch
                % backup variant
                disp('EEGLAB importer failed; falling back...');
                args = hlp_struct2varargin(opts,'suppress',{'samplerange'}, 'rewrite',{'timerange','blockrange'});
                res = pop_biosig(filename,args{:});
                opts.timerange = [];
                opts.channels = [];
            end
        case '.gdf'
            % General Data Format (.gdf)
            args = hlp_struct2varargin(opts,'suppress',{'samplerange'}, 'rewrite',{'timerange','blockrange'});
            res = pop_biosig(filename,args{:});
            opts.timerange = [];
            opts.channels = [];
        case '.mat'
            % BrainVision Analyzer Matlab file or BCI competition file
            try
                evalc('res = pop_loadbva(filename)');
                disp('Imported .mat file as a BrainVision Analyzer file.');
            catch
                % backup
                res = load(filename);
                % if this contains a single variable take that as the output
                if length(fieldnames(res)) == 1
                    res = struct2cell(res);
                    res = res{1};
                end
                % check if this is a BCI competition file
                if all(isfield(res,{'cnt','nfo'}))
                    disp('Parsing .mat file as a BCI competition MATLAB file.');
                    if isfield(res,'mrk')
                        % assemble a proper data set from the loaded pieces
                        evtypes = cellfun(@num2str,num2cell(res.mrk.y),'UniformOutput',false);
                        res = exp_eval(set_new('data',single(res.cnt'),'srate',res.nfo.fs,'event',struct('type',evtypes,'latency',num2cell(res.mrk.pos)),'chanlocs',res.nfo.clab));
                    else
                        res = exp_eval(set_new('data',single(res.cnt'),'srate',res.nfo.fs,'chanlocs',res.nfo.clab));
                    end
                end
            end
        case '.ds'
            % CTF folder
            optsctf = hlp_varargin2struct(opts,'trials',[]);
            res = pop_ctf_read(filename,opts.channels,opts.timerange,optsctf.trials);
            opts.channels = [];
            opts.timerange = [];
        case '.rdf'
            % ERPSS data
            res = pop_read_erpss(filename);
        case '.asc'
            % INStep ASC
            res = pop_loadascinstep(filename);
        case '.m4d'
            % 4D pdf file
            res = pop_read4d(filename);
        case '.dat'
            % BCI2000 file (by default with all runs concatenated)
            optsdat = hlp_varargin2struct(opts,{'mergeposition','MergePosition','merge_position'},true, ...
                {'concatruns','ConcatRuns','concat_runs'},true, {'maxevents','MaxEvents','max_events'},3000);
            res = BCI2000import(filename,false,optsdat.mergeposition,optsdat.concatruns,optsdat.maxevents);
        case '.xdf'
            % XDF files
            args = hlp_struct2varargin(opts);
            res = eeg_load_xdf(filename,args{:});
        case '.sto'
            % Storage files (specific to BCILAB)
            res = getfield(io_load(filename),'EEG');
        otherwise
            error('This file format has no known handler in BCILAB.');
    end
catch specific_error
    % the specific importers failed, fall back to the generic ones
    try
        % try to use FieldTrip
        warning off FieldTrip:unknown_filetype
        hdr = ft_read_header(filename);
        res = struct('setname','','filename','','filepath','','subject','','group','','condition','','session',[],'comments','','nbchan',hdr.nChans,...
            'trials',hdr.nTrials,'pnts',hdr.nSamples,'srate',hdr.Fs,'xmin',-hdr.nSamplesPre/hdr.Fs,'xmax',0,'times',[],'data',[],'icaact',[],'icawinv',[], ...
            'icasphere',[],'icaweights',[],'icachansind',[],'chanlocs',struct('labels',hdr.label),'urchanlocs',[],'chaninfo',[],'ref',[],'event',[],'urevent',[], ...
            'eventdescription',{{}}, 'epoch',[],'epochdescription',{{}},'reject',[],'stats',[],'specdata',[],'specicaact',[],'splinefile','','icasplinefile','', ...
            'dipfit',[],'history','','saved','no','etc',[]);
        if res.trials > 1
            error('This importer does not support epoched data.'); end
        if ~isempty(opts.timerange)
            opts.samplerange = 1+round((opts.timerange-res.xmin)*res.srate); end
        res.data = ft_read_data(filename,'chanindx',opts.channels,'begsample',opts.samplerange(1:end-1),'endsample',opts.samplerange(2:end));
        evt = ft_read_event(filename);
        res.event = struct('type',{evt.value},'category',{evt.type},'latency',{evt.sample},'duration',{evt.duration});
        res.xmax = res.xmin + (res.pnts-1)*res.srate;
        [opts.channels,opts.samplerange,opts.timerange] = deal([]);
    catch fileio_error
        try
            % try to use BioSig
            args = hlp_struct2varargin(opts,'suppress',{'channels','samplerange','timerange'});
            res = pop_biosig(filename,args{:});
        catch biosig_error
            disp('All possibly applicable importers for this file format failed.');
            disp('A list of error messages from the respectively tried loaders follows:');
            disp('    Error report of the format-specific loader:');
            env_handleerror(specific_error,6);
            disp('    Error report of the generic FileIO loader:');
            env_handleerror(fileio_error,6);
            disp('    Error report of the generic BioSig loader:');
            env_handleerror(biosig_error,6);
            fprintf('\n');
            error('BCILAB:io_loadset:cannot_load','Cannot load your data; please check for additional loader plugins at www.sccn.ucsd.edu/eeglab/plugins/.');
        end
    end
end

if allopts.casttodouble
    res.data = double(res.data);end

% if no markers present, automatically infer them
if isempty(res.event) || allopts.markerchannel.force_processing
    res = set_infer_markers('signal',res,allopts.markerchannel); end

% infer chanlocs fields (e.g. coordinates)
if allopts.infer_chanlocs
    res = set_infer_chanlocs(res,allopts.montage_disambiguation); end

% convert numeric event types to string
if isfield(res.event,'type') 
    numeric_mask = cellfun(@isnumeric,{res.event.type});
    if any(numeric_mask)
        disp('Converting all numeric event types to strings.');
        for k=find(numeric_mask)
            res.event(k).type = num2str(res.event(k).type); end
    end
end

% retain only channels with the selected type
if ~isempty(allopts.types)    
    if ischar(allopts.types)
        allopts.types = {allopts.types}; end
    matches = false;
    for t=1:length(allopts.types)
        matches = matches | strcmp({res.chanlocs.type},allopts.types{t}); end
    keep = find(matches);
    res.data = res.data(keep,:,:,:,:,:,:,:);
    res.chanlocs = res.chanlocs(keep);
    res.nbchan = size(signal.data,1);
end

% add data set meta-data
res.setname = allopts.setname;
res.filename = [name ext];
res.filepath = base;
res.subject = allopts.subject;
res.group = allopts.group;
res.condition = allopts.condition;
res.comments = allopts.comments;

% and reduce the data post-hoc if not already done so in the loader
if ~isempty(opts.channels)
    res = hlp_scope({'disable_expressions',true},@flt_selchans,res,{res.chanlocs(opts.channels).labels}); end
if ~isempty(opts.samplerange)
    res = hlp_scope({'disable_expressions',true},@set_selinterval,res,opts.samplerange,'samples'); end
if ~isempty(opts.timerange)
    res = hlp_scope({'disable_expressions',true},@set_selinterval,res,opts.timerange,'seconds'); end
if ~isempty(allopts.subsampled)
    res = hlp_scope({'disable_expressions',true},@flt_resample,res,allopts.subsampled);
    res.data = res.data(:,1+mod(((0:res.pnts-1) + round(res.etc.filter_delay*res.srate)),res.pnts-1));
end

if length(unique({res.chanlocs.labels})) ~= length(res.chanlocs)
    warning('bcilab:io_loadset:duplicate_channels','Multiple of your channels have the same label; this will likely give you errors during processing.'); end

% for faster processing further into the pipeline
res.tracking.timeseries_fields = {};

% add tracking information (note: the @rawdata expression indicates that this stage produces raw data...)
if isfield(res,'tracking') && isfield(res.tracking,'online_expression')
    % if the data set has already an online expression (e.g. processed data from a .set file), use that
    exp_endfun('set_online',res.tracking.online_expression)
else
    % otherwise we treat this stage as producing raw data
    exp_endfun('set_online',struct('head',@rawdata,'parts',{{{res.chanlocs.labels},unique({res.chanlocs.type})}}));
end
