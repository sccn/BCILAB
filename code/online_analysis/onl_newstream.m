function id = onl_newstream(name,varargin)
% Create a new data stream, and set up meta-data.
% Id = onl_newstream(Name, Options...)
%
% After a stream has been created, blocks of data can be appended to it using onl_append().
% Predictors can be linked to the stream using onl_newpredictor(), and their predictions (given the
% stream's most recent contents) can be queried using onl_predict(). This stream processing API is
% the basis for all online processing plugins to real-time experimentation environments (see
% code/online_plugins/*), as well as for the pseudo-online simulation function onl_simulate().
%
% A data stream is essentially a buffer of data (usually EEG, but may instead contain other sampled
% data, such as gaze coordinates, or motion capture coordinates). It has the same meta-data as a
% regular EEGLAB dataset, containing fields such as .srate, and .chanlocs, as well as special
% additional fields .smax, .buffer, and possibly others. All name-value pairs that are specified in
% the Options are added as fields to this data structure. Instead of .data, it has .buffer, which is
% a circular buffer that holds a segment of data ending with the most recently supplied sample.
%
% In:
%   Name : name of the stream; a variable of this name will be created in the workspace to
%          hold the stream's data.
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
%                * 'timestamps_len' : number of measures that should be averaged to yield a signal 
%                                     lag estimate, if time stamps are supplied online during onl_append
%                                     (default: 25)
%
%                * 'buffer_len' : maximum length of the signal buffer, in seconds (default: 10)
%                                 data will be lost if a predictor is not queried for longer than
%                                 this period                              
%
% Out:
%   Id : a unique id number for the predictor; same as name.predictorid
%
%
% Notes:
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
    error('The name of the online stream must be a valid variable name.'); end

% create the stream and add/reset user-specified fields
stream = exp_eval(set_new('buffer_len',10,'timestamps_len',25,'types',[],varargin{:}, ...
    'event',[],'urevent',[],'epoch',[],'icaact',[],'tracking',[]));

% check validity of user-specified fields
if isempty(stream.srate)
    error('A sampling rate must be specified via the ''srate'' option.'); end
if isempty(stream.chanlocs)
    error('Channel info must be specified via the ''chanlocs'' option.'); end
if isempty(stream.types)
    stream.types = unique({stream.chanlocs.type}); end
stream.nbchan = length(stream.chanlocs);

% add buffer fields
stream.buffer_len = stream.buffer_len*stream.srate;     % max. number of samples in the stream's buffer (the buffer holds a running sub-range of the entire stream sample range)
stream.buffer = zeros(stream.nbchan,stream.buffer_len); % circular data buffer
stream.smax = 0;                                        % index of last sample in the stream; 
                                                        % smin = 1+max(0,smax-buffer_len)
                                                        % # of valid samples in buffer = smax-smin+1
                                                        % current data range:
                                                        % buffer(:,1+mod((smin:smax)-1,buffer_len))

stream.timestamps = zeros(stream.timestamps_len,2); % circular array of time stamp estimates
stream.timestamps_ptr = 0;                          % last write pointer into timestamps

% remove fields that we do not track
stream = rmfield(stream,{'pnts','xmax','data'});
if stream.xmin ~= 0
    % put xmin into timestamps so that it blends with any later time stamp updates during onl_append
    stream.timestamps_ptr = 1+mod(stream.timestamps_ptr,stream.timestamps_len);
    stream.timestamps(stream.timestamps_ptr,:) = [stream.xmin 0];
end

% assign a unique id to the stream (for users of the stream)
id = fresh_id('bcilab_streams');
stream.streamid = id;

% add a few different time stamps
try
    stream.xmin_micro = tic;
    stream.xmin_mili = java.lang.System.currentTimeMillis();
    stream.xmin_nano = java.lang.System.nanoTime();
catch, end

% create a new stream and place it in the base workspace
assignin('base',name,stream);
