function run_writedataset(varargin)
% Output a raw online stream into an EEGLAB dataset.
% run__writedataset(SourceStream,FileName,UpdateFrequency,StartDelay)
%
% This function does not do any processing, but just saves a stream to a file (possibly in parallel
% to some other operation processing it).
%
% In:
%   SourceStreamName : Optional name of the stream data structure in the MATLAB base workspace to
%                      take as the data source (previously created with onl_newstream).
%                      (default: 'laststream')
%
%   FileName : File name to write to (default: 'lastdata.set')
%
%   UpdateFrequency : The are at which new chunks of data are appended to the file, in Hz (default: 1)
%
%   StartDelay : Start-up delay before real-time operation begins; grace period until file is being
%                written to, in s. (default: 3)
%
% Examples:
%   % write an input stream (named 'mystream') to a file named 'recording.set' (EEGLAB dataset)
%   run_writedataset('mystream','recording.set')
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

declare_properties('name','File');

% define arguments
arg_define(varargin, ...
    arg({'in_stream','SourceStreamNames','SourceStream'}, 'laststream',[],'Input Matlab stream name(s). Optional names of stream data structures in the MATLAB base workspace to consider as possible data sources (previously created with onl_newstream); if a stream contains all channels that are needed by the predictor, or alternatively has the right number and type of channels it will be considered as a potential source stream unless ambiguous.'), ...
    arg({'out_filename','FileName'},'lastdata.set',[],'The file name to write to.'), ...
    arg({'update_freq','UpdateFrequency'},1,[0 Inf],'Update frequency. This is the rate at which new chunks of data are appended to the file.'), ...
    arg({'start_delay','StartDelay'}, 3, [0 Inf],'Start-up delay. Delay before real-time operation begins; grace period until file is written.'));

out_filename = env_translatepath(out_filename);

% open the stream and write the initial set file header...
stream = evalin('base',in_stream);
% create missing fields
stream.data = randn(stream.nbchan,1024);
stream.pnts = size(stream.data,2);
stream.xmax = stream.xmin + (stream.pnts-1)/stream.srate;
[fp fn fe] = fileparts(out_filename);
% remove superfluous fields
eeg = rmfield(stream,{'buffer','smax','buffer_len','timestamps','timestamps_len','timestamps_ptr','streamid'});
stream.timestamp_at_beginning = toc(uint64(0));
eeg = pop_saveset(eeg,'filename',[fn fe],'filepath',quickif(isempty(fp),env_translatepath('bcilab:/userdata'),fp),'savemode','twofiles');
% re-create the fdt file...
delete(fullfile(eeg.filepath, eeg.datfile));
fid = fopen(fullfile(eeg.filepath, eeg.datfile),'wb','ieee-le');
if fid == -1
    error('Cannot write output file, check permission and space.'); end;
% create timer (which periodically writes to the stream)
t = timer('ExecutionMode','fixedRate', 'Name',[in_stream '_write_timer'], 'Period',1/update_freq, ...
    'StartDelay',start_delay, 'TimerFcn',@(obj,varargin) append_data(in_stream,fid,stream.streamid,obj,eeg));

% start timer
start(t);


% timer callback: visualization
function append_data(stream,fid,streamid,timerhandle,eeg)
try
    % check if the stream and the predictor are still there
    s = evalin('base',stream);
    if s.streamid ~= streamid
        error('Stream changed.'); end

    % get an updated chunk of data
    samples_to_get = min(s.buffer_len, s.smax-ftell(fid)/(4*s.nbchan));
    chunk = s.buffer(:, 1+mod(s.smax-samples_to_get:s.smax-1,s.buffer_len));    
    
    % and write it into the file
    fwrite(fid,chunk,'float');    
catch e
     if ~strcmp(e.identifier,'MATLAB:UndefinedFunction')
        env_handleerror(e); end
    finalize_dataset(fid,eeg);
    % interrupted: make sure that the file gets closed
    stop(timerhandle);
    delete(timerhandle);
end


function finalize_dataset(fid,EEG)
samples = ftell(fid)/(4*EEG.nbchan);
fclose(fid);
EEG.pnts = samples;
EEG.data = EEG.datfile;
EEG.xmax = EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.timestamp_at_end = toc(uint64(0));

save(fullfile(EEG.filepath, EEG.filename), '-v6', '-mat', 'EEG');
