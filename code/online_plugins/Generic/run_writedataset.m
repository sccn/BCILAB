function run_writedataset(varargin)
% Output a raw online stream into an EEGLAB dataset.
% run__writedataset(SourceStream,FileName,UpdateFrequency,StartDelay)
%
% This function does not do any processing, but just saves a stream to a file (possibly in parallel
% to some other operation processing it).
%
% In:
%   SourceStream : real-time stream name to read from (in MATLAB workspace) (default: 'laststream')
%
%   FileName : File name to write to (default: 'lastdata.set')
%
%   UpdateFrequency : update frequency, in Hz (default: 1)
%
%   StartDelay : Start-up delay before real-time processing begins; grace period until file is being
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
    arg({'in_stream','SourceStream'}, 'laststream',[],'Input online stream. This is the stream that shall be written to disk.'), ...
    arg({'out_filename','FileName'},'lastdata.set',[],'The file name to write to.'), ...
    arg({'update_freq','UpdateFrequency'},1,[],'Update frequency. This is the rate at which data is written.'), ...
    arg({'start_delay','StartDelay'}, 3, [],'Start-up delay. Delay before real-time processing begins; grace period until file is written.'));

out_filename = env_translatepath(out_filename);

% open the stream and write the initial set file header...
stream = evalin('base',in_stream);
% create missing fields
stream.data = randn(stream.nbchan,1024);
stream.pnts = size(stream.data,2);
stream.xmax = stream.xmin + (stream.pnts-1)/stream.srate;
[fp fn fe] = fileparts(out_filename);
% remove superfluous fields
eeg = rmfield(stream,{'buffer','smax','buffer_len','timestamps','timestamps_len','tmax','streamid'});
stream.timestamp_at_beginning = toc(uint64(0));
eeg = pop_saveset(eeg,'filename',[fn fe],'filepath',fastif(isempty(fp),env_translatepath('bcilab:/userdata'),fp),'savemode','twofiles');
% re-create the fdt file...
delete(fullfile(eeg.filepath, eeg.datfile));
fid = fopen(fullfile(eeg.filepath, eeg.datfile),'wb','ieee-le');
if fid == -1
    error('Cannot write output file, check permission and space.'); end;

% create a temporary text file to store stream event markers
event_filename = env_translatepath(['bcilab:/userdata/', fn, '_events.txt']);
marker_fid = fopen(event_filename, 'w+');
if marker_fid == -1
    error('Cannot write marker output file, check permission and space.');
end
% write a basic header to the event marker file
fprintf(marker_fid, '%s\t%s\n', 'type', 'latency');

% create timer (which periodically writes to the stream)
t = timer('ExecutionMode','fixedRate', 'Name',[in_stream '_write_timer'], 'Period',1/update_freq, ...
    'StartDelay',start_delay, 'TimerFcn',@(obj,varargin) append_data(in_stream,fid,marker_fid,stream.streamid,obj,eeg));

% start timer
start(t);


% timer callback: visualization
function append_data(stream,fid,marker_fid,streamid,timerhandle,eeg)
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
    
    % check if this data chunk contains any event markers
    marker_chunk = s.marker_pos(:, 1+mod(s.smax-samples_to_get:s.smax-1,s.buffer_len));
        
    % if so, write them into the marker file
    if nnz(marker_chunk) > 0
        % find the sample(s) in this chunk with events
        [marker_pos_in_sample, marker_pos_in_chunk] = find(marker_chunk);

        % EEG samples which correspond to positions in the chunk
        chunk_samples = s.smax - (samples_to_get - 1) : s.smax;
        
        for m = 1:length(marker_pos_in_chunk)
            marker_idx = marker_chunk(marker_pos_in_sample(m), marker_pos_in_chunk(m));
            marker_type = s.marker_buffer(marker_idx).type;
            marker_sample_whole = chunk_samples(marker_pos_in_chunk(m));
            marker_sample_fractional = s.marker_buffer(marker_idx).latency;
            marker_sample = marker_sample_whole + marker_sample_fractional;
            fprintf(marker_fid, '%s\t%i\n', marker_type, marker_sample);
        end 
    end
    
catch e
     if ~strcmp(e.identifier,'MATLAB:UndefinedFunction')
        env_handleerror(e); end
    finalize_dataset(fid,marker_fid,eeg);
    % interrupted: make sure that the file gets closed
    stop(timerhandle);
    delete(timerhandle);
end


function finalize_dataset(fid,marker_fid,EEG)
samples = ftell(fid)/(4*EEG.nbchan);
fclose(fid);
EEG.pnts = samples;
EEG.data = EEG.datfile;
EEG.xmax = EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.timestamp_at_end = toc(uint64(0));

% load the event marker file from disk and include the markers in the EEG .set
marker_file = [EEG.filepath, filesep, regexprep(EEG.filename, '.set', '_events.txt')];
frewind(marker_fid);
marker_data = textscan(marker_fid, '%s%d', 'Delimiter', '\t', 'HeaderLines', 1);
if ~isempty(marker_data{1})
   EEG.event = struct('type', marker_data{1}', 'latency', num2cell(marker_data{2}'), ...
        'urevent', num2cell((1:length(marker_data{1}))));
end
fclose(marker_fid);
delete(marker_file);

save(fullfile(EEG.filepath, EEG.filename), '-v6', '-mat', 'EEG');



