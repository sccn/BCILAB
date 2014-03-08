function id = onl_newstream(name,varargin)
% Create a new data stream, and set up meta-data.
% Id = onl_newstream(StreamName, Options...)
%
% After a stream has been created, blocks of data can be appended to it using onl_append().
% Predictors can be linked to the stream using onl_newpredictor(), and their predictions (given the
% stream's most recent contents) can be queried using onl_predict(). One can alternatively also get
% the most recent k seconds of raw data via onl_peek() or set up a filter pipeline using
% onl_newpipeline() and then query the most recent k samples of filtered data using onl_filtered().
% This stream API is the basis for all online processing plugins (see code/online_plugins/*), as
% well as for the pseudo-online simulation function onl_simulate().
%
% A stream has almost the same layout as an EEGLAB dataset struct; particularly, all the meta-data
% (such as .srate and .chanlocs) is exactly the same. The difference is that a) the stream holds
% only the last n seconds of data and markers that were appended (since all read operations on it
% refer to portions of the most recent data), and b) markers and data are not in the regular EEGLAB
% format for efficiency reasons (see Notes for internal details).
%
% All name-value pairs that are specified in the Options are added as fields to this data structure.
%
% In:
%   StreamName : name of the stream; a variable with this name will be created in the MATLAB workspace 
%                to hold the stream's data. If such a variable already exists it will be overridden.
%
%   Options... : meta-data fields that are to be included in the stream, given as name-value pairs. 
%                A struct with meta-data (e.g., a dataset), may also be specified.
%
%                * 'srate' : the sampling rate, in Hz, of the stream (mandatory)
%
%                * 'chanlocs' : per-channel meta-data ('locs' is due to the EEGLAB heritage);
%                               struct array with fields
%                                 .labels : channel name (mandatory)
%                                 .type   : channel type (typical types are EEG, EMG, EOG, ECG, GSR,
%                                           Gaze, ECoG, NIRS, MoCap)
%                                 additional domain-specific fields (e.g., X,Y,Z coordinates for EEG 
%                                 channels may be added)
%                               can also be a cell array of strings, which are then assumed to 
%                               represent the channel labels
%
%                * 'xmin' : optional time, in seconds, of the first sample in some arbitrary clock 
%                           domain; if non-zero xmin's are specified for multiple streams that are 
%                           being processed by the same predictor, they must be all in the same
%                           domain. If time stamps are supplied during onl_append the xmin/xmax 
%                           are re-adapted consistently with the time stamps. (default: 0)
%
%                * 'types' : cell-string array of type(s) of the stream data (default: union of the 
%                            channel types); this can be utilized to declare that a stream carries
%                            data of some types even though no proper channel labels are available,
%                            which allows for lax fallback matching of streams and filter chains
%                            of predictors in the absence of channel labels
%
%                * domain-specific additional fields, such as icawinv, dipfit, or ref
%
%                Additional online-specific fields:
%
%                * 'buffer_len' : maximum length of the signal buffer, in seconds (default: 10).
%
%                * 'marker_buffer_len' : maximum number of marker records held in the marker buffer
%                                       (default: 1000).
%
%                * 'timestamps_len' : number of time stamps that should be averaged to yield a signal 
%                                     lag estimate, if time stamps are supplied online during onl_append
%                                     (default: 25)
%
% Out:
%   Id : a unique id number for the predictor; same as StreamName.streamid
%
%
% Notes:
%   For efficiency, the storage layout for the stream's data is not an array that goes from the
%   first sample to the last, but rather a circular buffer where data wraps around; to avoid
%   confusion, the stream has empty .data and .event fields. Instead of .data, it has a .buffer
%   field (a circular buffer that holds the most recently added n seconds of data) and .smax (the
%   number of samples appended so far, which is used as an index into the ring buffer, using
%   wrap-around indexing). Instead of .event it has a field .marker_buffer (circular buffer of
%   marker structs), .mmax (number of marker records appended so far and serving as wrap-around
%   index into .marker_buffer), and lastly .marker_pos (sparse time series that holds indices into
%   .marker_buffer; like .buffer, the .marker_pos time series is a circular buffer). Please avoid
%   using these fields directly -- instead always use the functions onl_append and onl_peek to
%   read/write to the stream since the internal format may change (also, fields need to be accessed
%   in the correct order to support concurrent read/write operations from timers).
%
%   The stream can be distinguished from a regular EEGLAB dataset by having the additional fields
%   'buffer' and 'smax'. smax indicates how many samples have been appended to the stream so far 
%   (though only the most recent subset of them is held in the buffer).
%
%   If 'clear' is run in the workspace, the stream will be erased.
%
% Example:
%   % calibrate (or load) a model
%   mydata = io_loadset('data:/projects/attention/mason.bdf');
%   [loss,model] = bci_train({'data',mydata, 'approach',myapproach});
%   
%   % create a new online stream
%   onl_newstream('mystream','srate',200,'chanlocs',{'C3','Cz','C4'});
%   % create a new predictor from a model
%   onl_newpredictor('mypredictor',model);
%
%   while 1
%       % obtain a new block of data from some acquisition system
%       datachunk = get_new_samples();
%       % feed the block into the online stream
%       onl_append('mystream',datachunk);
%       % obtain predictions
%       outputs = onl_predict('mypredictor');
%       % display them
%       bar(outputs); drawnow; pause(0.1);
%   end
%
% See also:
%   onl_newpredictor, onl_append, onl_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% check validity of name
if ~exist('name','var')
    error('Please specify a name for the stream.'); end
if ~isvarname(name)
    error('The given StreamName argument must be a valid variable name, but was: %s',hlp_tostring(name,10000)); end

% create the stream and add/reset user-specified fields
stream = exp_eval(set_new('buffer_len',10,'timestamps_len',25,'marker_buffer_len',1000,'types',[],varargin{:}, ...
    'event',[],'urevent',[],'epoch',[],'icaact',[],'tracking',[]));

% check validity of user-specified fields
if isempty(stream.srate)
    error('A sampling rate must be specified via the ''srate'' option.'); end
if ~isnumeric(stream.srate) || ~isscalar(stream.srate) || stream.srate == 0
    error('The given sampling rate is invalid: %s.',hlp_tostring(stream.srate,10000)); end
if isempty(stream.chanlocs)
    error('Channel info must be specified via the ''chanlocs'' option.'); end
if ~isfield(stream.chanlocs,'labels')
    error('The given chanlocs field is lacking the required .labels field.'); end
if isempty(stream.types)
    stream.types = unique({stream.chanlocs.type}); 
elseif ~iscellstr(stream.types)
    error('The given types argument must be a cell array of strings, but was: %s',hlp_tostring(stream.types,10000));
end
stream.nbchan = length(stream.chanlocs);
if stream.buffer_len > 120
    disp_once('Warning: The given stream buffer length is very long and might cause inefficiencies in processing.'); end
if stream.buffer_len < 1
    disp_once('Warning: The given stream buffer length is very short (typically at least a few seconds).'); end
if stream.marker_buffer_len > 100*stream.buffer_len
    disp_once('Warning: The given marker buffer length is very long and might cause inefficiencies in processing (100 per each second in the data buffer is a reasonable upper limit).'); end
if stream.marker_buffer_len < 3*stream.buffer_len
    disp_once('Warning: The given marker buffer length is rather short (typically at least 3 per each second in the data buffer).'); end
if stream.timestamps_len > 100
    disp_once('Warning: The given timestamps_len is very long and might cause inefficiencies in processing (100 is a good upper limit).'); end
if stream.timestamps_len < 10
    disp_once('Warning: The given timestamps_len is rather short and might cause jitter in the time-stamp estimates (10 is a reasonable lower limit).'); end

% add data buffer fields
stream.buffer_len = stream.buffer_len*stream.srate;     % max. number of samples in the stream's buffer (the buffer holds a running sub-range of the entire stream sample range)
stream.buffer = zeros(stream.nbchan,stream.buffer_len); % circular data buffer
stream.smax = 0;                                        % index of last sample in the stream; 
                                                        % smin = 1+max(0,smax-buffer_len)
                                                        % # of valid samples in buffer = smax-smin+1
                                                        % current data range:
                                                        % buffer(:,1+mod((smin:smax)-1,buffer_len))

% add marker buffer fields                                                        
max_simultaneous_markers = 1000; % note: this constant must equal the one in onl_append
stream.marker_pos = sparse(zeros(max_simultaneous_markers,stream.buffer_len));  % sparse array mapping from time (sample indices) to marker_buffer; indexed by [#sample,#duplicate]
stream.marker_buffer = struct(...                                               % buffer of past k marker records; latency is the fractional offset within the sample (not absolute)
    'type',cell(1,stream.marker_buffer_len), ...      
    'latency',cell(1,stream.marker_buffer_len));                               
stream.mmax = 0;

% add time stamp fieds
stream.timestamps = zeros(stream.timestamps_len,2); % circular array of time stamp estimates
stream.tmax = 0;                                    % number of timestamps appended so far

% remove fields that we do not track
stream = rmfield(stream,{'pnts','xmax','data'});
if stream.xmin ~= 0
    % put xmin into timestamps so that it blends with any later time stamp updates during onl_append
    stream.tmax = 1+mod(stream.tmax,stream.timestamps_len);
    stream.timestamps(stream.tmax,:) = [stream.xmin 0];
end

% assign a unique id to the stream (for users of the stream)
id = fresh_id('bcilab_streams');
stream.streamid = id;

% add a few different time stamps
try
    stream.xmin_micro = tic;
    stream.xmin_mili = java.lang.System.currentTimeMillis();
    stream.xmin_nano = java.lang.System.nanoTime();
catch
end

% create a new stream and place it in the base workspace
assignin('base',name,stream);
