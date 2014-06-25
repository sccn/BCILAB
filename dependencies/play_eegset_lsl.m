function handle = play_eegset_lsl(dataset,datastreamname,eventstreamname,looping,background,update_interval)
% Play back a continuous EEGLAB dataset over LSL.
% Handle = play_eegset_lsl(Dataset,DataStreamName,EventStreamName,Looping,Background,UpdateInterval)
%
% In:
%   Dataset : EEGLAB dataset struct to play.
%
%   DataStreamName : Name of the data stream to create (default: 'EEGLAB')
%
%   EventStreamName : Name of the event stream to create (default: 'EEGLAB-Markers')
%
%   Looping : Whether play back the data in a loop (default: true)
%
%   Background : Whether to run in the background; see example. (default: false)
%
%   UpdateInterval : Interval between updates, in s (default: 0)
%
% Out:
%   Handle : A handle that can be used to stop background playback, by calling stop(myhandle).
%
% Examples:
%   EEG = io_loadset('data:/tutorial/flanker_task/12-08-001_ERN.vhdr');
%
%   % Play a recording (blocking, press Ctrl+C to cancel)
%   play_eegset_lsl(EEG);
%
%   % Play a recording in the background
%   h = play_eegset_lsl(EEG,[],[],[],true);
%   % stop playback after some time
%   stop(h);
%
%   % Play with looping disabled
%   play_eegset_lsl(EEG,[],[],false);
%
%   Play with a different name for the data stream and marker stream
%   play_eegset_lsl(EEG,'MyDataStream','MyMarkerStream');
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-23

    if ~exist('datastreamname','var') || isempty(datastreamname)
        datastreamname = 'EEGLAB'; end
    if ~exist('eventstreamname','var') || isempty(eventstreamname)
        eventstreamname = 'EEGLAB-Markers'; end
    if ~exist('looping','var') || isempty(looping)
        looping = true; end
    if ~exist('background','var') || isempty(background)
        background = false; end
    if ~exist('update_interval','var') || isempty(update_interval)
        update_interval = 0; end
    

    if background && ~nargout
        error('To play back in the background you need to specify a return value (otherwise you playback could not be stopped).'); end
    
    % evaluate the dataset if necessary
    if ~isfield(dataset,'data') && all(isfield(dataset,{'head','parts'})) && exist('exp_eval','file')
        dataset = exp_eval(dataset); end
    
    disp('Loading LSL library...');
    lib = lsl_loadlib();

    disp('Creating data streaminfo...');
    datainfo = lsl_streaminfo(lib,datastreamname,'EEG',dataset.nbchan,dataset.srate,'cf_float32',matrix_hash(dataset.data));
    desc = datainfo.desc();
    channels = desc.append_child('channels');
    for c=1:length(dataset.chanlocs)
        channel = channels.append_child('channel');
        channel.append_child_value('label',dataset.chanlocs(c).labels);
        if isfield(dataset.chanlocs,'type')
            channel.append_child_value('type',dataset.chanlocs(c).type); end
        if all(isfield(dataset.chanlocs,{'X','Y','Z'}))
            location = channel.append_child('location');
            if isfield(dataset,'chaninfo') && isfield(dataset.chaninfo,'nosedir') && strcmp(dataset.chaninfo.nosedir,'+X')
                location.append_child_value('X',dataset.chanlocs(c).Y*1000);
                location.append_child_value('Y',dataset.chanlocs(c).X*1000);
                location.append_child_value('Z',dataset.chanlocs(c).Z*1000);
            else
                location.append_child_value('X',dataset.chanlocs(c).X*1000);
                location.append_child_value('Y',dataset.chanlocs(c).Y*1000);
                location.append_child_value('Z',dataset.chanlocs(c).Z*1000);
            end
        end
    end

    disp('Opening data outlet...');
    dataoutlet = lsl_outlet(datainfo);

    disp('Creating marker streaminfo...');
    markerinfo = lsl_streaminfo(lib,eventstreamname,'Markers',1,0,'cf_string',[matrix_hash(dataset.data) '-markers']);
    disp('Opening marker outlet...');
    markeroutlet = lsl_outlet(markerinfo);

    % make a sparse index map for faster event lookup
    disp('Preparing marker map...');
    if ~isempty(dataset.event)
        marker_types = {dataset.event.type};
        [marker_map,residuals] = sparse_binning(min(dataset.pnts,max(1,[dataset.event.latency])),[],dataset.pnts);
    else
        marker_map = sparse(1,dataset.pnts);
    end

    last_pos = 0;
    start_time = tic;
    disp('Now playing back...');
    if background
        handle = timer('ExecutionMode','fixedRate', 'Period',0.01, 'TasksToExecute',fastif(looping,Inf,dataset.xmax/0.01), 'TimerFcn',@update, 'StopFcn',@(t,varargin)delete(t));
        start(handle);
    else
        while looping || last_pos<dataset.pnts
            if ~update()
                pause(0.001); end
            if update_interval > 0
                pause(update_interval); end
        end
        disp('Reached end of data set.');
    end


    function need_update = update(varargin)
        % perform an update (send new data over LSL, if any)
        pos = round(1+toc(start_time)*dataset.srate);
        need_update = pos > last_pos;
        if need_update
            range = (last_pos+1):pos;
            wraprange = 1+mod(range-1,dataset.pnts);
            now = lsl_local_clock(lib);
            % push markers in range
            [ranks,sample_indices,marker_indices] = find(marker_map(:,wraprange));
            if any(ranks)
                % calculate time stamps
                marker_times = now + (sample_indices(:) + vec(residuals(marker_indices)) - length(range))/dataset.srate;
                % send them off
                for m=1:length(marker_times)
                    markeroutlet.push_sample(marker_types(marker_indices(m)),marker_times(m)); end
            end
            % push data chunk
            dataoutlet.push_chunk(double(dataset.data(:,wraprange)),now);
            last_pos = pos;
        end
    end

    function stop_playback(obj,varargin)
        delete(obj);
        disp('Playback stopped.');
    end

end

    
function hash = matrix_hash(data)
    % Calculate a hash of a given matrix.
    max_java_memory = 2^26; % 64 MB
    hasher = java.security.MessageDigest.getInstance('MD5');
    data = typecast(full(data(:)),'uint8');
    if length(data) <= max_java_memory
        hasher.update(data);
    else
        numsplits = ceil(length(data)/max_java_memory);
        for i=0:numsplits-1
            range = 1+floor(i*length(data)/numsplits) : min(length(data),floor((i+1)*length(data)/numsplits));
            hasher.update(data(range));
        end
    end
    hash = dec2hex(typecast(hasher.digest,'uint8'),2);
    hash = char(['X' hash(:)']);
end


function [B,R] = sparse_binning(V,rows,columns)
    % Round and bin the given values into a sparse array.
    % [BinIndices,Residuals] = sparse_binning(Values,Padding)
    % 
    % In:
    %   Values : vector of values that yield valid indices (>1) when rounded.
    %
    %   Rows : number of rows to reserve (default: []=as many as needed)
    %
    %   Columns : number of columns to reserve (default: []=as many as needed)
    %
    % Out:
    %   BinIndices : sparse array containing indices into Values, where the horizontal axis is the
    %                rounded value and the vertical axis are the ranks of values that fall into the same
    %                bin.
    %
    %   Residuals : optionally the fractional offset between the values and the bin indices.

    [V,order] = sort(V(:)'); %#ok<UDIM,TRSRT>
    % get the bin index where each value falls
    bins = round(V);
    % get the within-bin rank of each value
    ranks = [false ~diff(bins)];
    ranks = 1+cumsum(ranks)-cummax(~ranks.*cumsum(ranks));
    % create a sparse matrix whose k'th column holds the indices to values in the k'th bin
    if nargin == 1
        B = sparse(ranks,bins,order);
    elseif nargin == 2
        B = sparse(ranks,bins,order,rows,max(bins));
    else
        if isempty(rows)
            rows = max(ranks); end
        if isempty(columns)
            columns = max(bins); end
        B = sparse(ranks,bins,order,rows,columns);
    end
    if nargout>0
        R = V-bins; end
end


function v = cummax(v)
    % Calculate cumulative maximum of a vector.
    m = v(1);
    for k = 2:numel(v)
        if v(k) <= m
            v(k) = m;
        else
            m = v(k);
        end
    end
end

function x = vec(x)
    % Vectorize a given array.
    x = x(:);
end
