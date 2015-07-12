function onl_append(name, chunk, markers, stamp)
% Append a block of raw data to a stream.
% onl_append(StreamName, Chunk, Markers, Timestamp)
%
% This function appends a chunk of data (#samples x #channels) and possibly some markers to an
% online stream that was previously created with onl_newstream. It is also possible to pass in the
% time stamp of the most recent sample for accurate multi-stream synchronization (but usually not
% necessary).
%
% The number of columns in the appended data block must match the number of channels in the stream,
% and the data must have been sampled at the rate specified in the stream's .srate field.
%
% The sequence of samples fed into a stream (over the course of multiple calls to onl_append) should
% have no omissions (i.e., samples should not be skipped if possible) or repetitions (i.e., samples
% that have been supplied in a previous call must not be fed again) since the filters would then see
% jumps and glitches in the data as a result.
%
% Processing (i.e., appending samples and calculating predictions from them) does not have to happen
% in real time (i.e., the actual wall-clock time does not affect the calulated results).
%
% In:
%   StreamName : Name of the stream data structure in the MATLAB workspace to which the data should 
%                be appended (must have been previously created with onl_newstream).
%
%   Chunk : [#Channels x #Samples] matrix of raw data
%
%   Markers : Optional struct array of markers in the chunk with one struct per marker. Each
%             struct must have a .type field that contains the marker type (string) and a .latency
%             field that contains the marker latency (where 1 is the first sample in the chunk and
%             #Samples is the last sample in the chunk). Additional fields are allowed, but note
%             that the field names of multiple successively appended marker structs must be
%             identical. Marker latencies should never exceed the bounds of the chunk. (default: empty)
%
%   Timestamp : Optional time stamp for the last sample in the chunk, in seconds, in some arbitary
%               time domain which is consistent across all streams used to make a given prediction.
%               (default: computed from xmin)
%
% Notes:
%   For many applications, accurate time codes are not necessary for reasonably good performance (or
%   not worth the development burden). For some applications, accurate timing across streams is
%   essential (especially if the relative lags are significant and temporal structure across streams
%   is being used for predictions). In this case, the time codes should be reasonably precise
%   estimates of the time stamp of the last sample, in a time domain that is the same across all
%   streams that are being used by a given predictor. This can be seconds since startup, or, for
%   example, the estimated (negative) age of the sample.
%
% Examples:
%   % get a chunk of new data from a device and append it to a (previously created) online stream
%   mychunk = get_data_from_device();
%   onl_append('mystream',mychunk);
%
%   % append both a 32-channel/1000-sample chunk and 3 markers
%   mychunk = randn(32,1000);
%   mymarkers = struct('type',{'test','sdfsdf','asd'},'latency',{10,300,801.5});
%   onl_append('mystream',mychunk,mymarkers);
%
% See also:
%   onl_newstream, onl_newpredictor, onl_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03

max_simultaneous_markers = 1000; % note: this constant must equal the one in onl_newstream

try
    if isempty(chunk)
        return; end
    
    % append markers
    if nargin>3 || (nargin==3 && isstruct(markers) && ~isempty(markers))        
        % create a sparse index array for the markers that can be streamed efficiently
        [marker_pos,residual] = sparse_binning(max(0.5,min(size(chunk,2)+0.5,[markers.latency])),max_simultaneous_markers,size(chunk,2));
        % rewrite the marker latency into fractional offset relative to the marker index
        [markers.latency] = arraydeal(residual);
        % append the marker chunk to the marker buffer
        perform_update(name,'marker_chunk',markers, ...
            ['try\n' ...
            '    X.marker_buffer(:,1+mod(X.mmax:X.mmax+size(X_marker_chunk,2)-1,X.marker_buffer_len)) = X_marker_chunk;\n' ...
            'catch e\n' ...
            '    if any(strcmp(e.identifier,{''MATLAB:heterogeneousStrucAssignment'',''MATLAB:heterogenousStrucAssignment''}))\n' ...
            '        X.marker_buffer = X_marker_chunk(ones(1,X.marker_buffer_len));\n' ...
            '        X.marker_buffer(:,1+mod(X.mmax:X.mmax+size(X_marker_chunk,2)-1,X.marker_buffer_len)) = X_marker_chunk;\n' ...
            '    else\n' ...
            '        rethrow(e);\n' ...
            '    end\n' ...
            'end']);
        % append the marker positions
        perform_update(name,'marker_pos',marker_pos, [...
            'X.marker_pos(:,1+mod(X.smax:X.smax+size(X_marker_pos,2)-1,X.buffer_len)) = X_marker_pos+logical(X_marker_pos)*X.mmax;\n'...
            'X.mmax = X.mmax + nnz(X_marker_pos);']);
        % append data chunk
        perform_update(name,'chunk',chunk,[...
            'X.buffer(:,1+mod(X.smax:X.smax+size(X_chunk,2)-1,X.buffer_len)) = X_chunk;\n' ...
            'X.smax = X.smax + size(X_chunk,2);']);
    else
        % append data chunk and clear marker positions in the range
        perform_update(name,'chunk_clr',chunk,[...
            'X_range = 1+mod(X.smax:X.smax+size(X_chunk_clr,2)-1,X.buffer_len);\n' ...
            'X.marker_pos(:,X_range) = 0;\n' ....
            'X.buffer(:,X_range) = X_chunk_clr;\n' ...
            'X.smax = X.smax + size(X_chunk_clr,2);']);
    end

    % update time stamps
    if (nargin==3 && isnumeric(markers)) || nargin>3
        if nargin == 3
            stamp = markers; 
            clear markers;
        end
        if isempty(stamp)
            return; end
        perform_update(name,'stamp',stamp,[...
            ... % insert timestamp at the end if the buffer (contains the time in timestamp domain and in
            ... % 0-based stream time domain)
            'X.tmax = 1+mod(X.tmax,X.timestamps_len-1);' ...
            'X.timestamps(X.tmax,:) = [X_stamp X.smax/X.srate];' ...
            ... % estimate time lag (offset of time stamp time relative to 0-based stream time)
            'X.lag = sum(X.timestamps(:,1) - X.timestamps(:,2))/nnz(X.timestamps(:,2));' ...
            ... % derive xmin/xmax in the time stamp domain
            '[X.xmin,X.xmax] = deal((1+max(0,X.smax-X.buffer_len))/X.srate + X.lag,X.smax/X.srate + X.lag);']);
    end
catch e
    % diagnose the error
    try        
        % validate name argument
        if nargin < 2
            error('At least the Name and Chunk argument must be given to onl_append.'); end
        if ~ischar(name) || ~isempty(name)
            error('The given Name argument must be a string that refers to a variable in the MATLAB workspace, but was: %s',hlp_tostring(name,10000)); end
        if ~isvarname(name)
            error('The given Name argument must be a valid variable name in the MATLAB workspace, but was: %s',name); end        
        try
            stream = evalin('base',name);
        catch %#ok<CTCH>            
            error('A stream named %s does not exist in the MATLAB workspace.',name);
        end

        % validate chunk argument
        if ndims(chunk) ~= 2 %#ok<ISMAT>
            error('The given Chunk argument must be a [#Channels x #Samples] matrix, but was: %s',hlp_tostring(chunk,10000)); end
        if ~isnumeric(chunk)
            error('The given Chunk argument must be a numeric matrix, but was: %s',hlp_tostring(chunk,10000)); end
        
        % validate stream fields
        if ~isstruct(stream) || ~isscalar(stream)
            error('The stream %s is not a 1x1 struct, but was: %s. It might have been overwritten by another function or you use the wrong stream name.',name,hlp_tostring(stream,10000)); end
        if ~isfield(stream,'buffer')
            error('The stream %s has no .buffer field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if ~isnumeric(stream.buffer)
            error('The stream %s has an invalid buffer. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if ~isfield(stream,'buffer_len')
            error('The stream %s has no buffer_len field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if ~isscalar(stream.buffer_len) || ~isnumeric(stream.buffer_len)
            error('The stream %s has an invalid buffer_len field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if size(stream.buffer,2) ~= stream.buffer_len
            error('The stream %s has a buffer with an invalid length. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if ~isfield(stream,'smax')
            error('The stream %s has no smax field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if ~isscalar(stream.smax) || ~isnumeric(stream.smax)
            error('The stream %s has an invalid smax field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        if size(stream.buffer,1) ~= size(chunk,1)
            error('Number of channels in the supplied chunk (%i) does not match the number of channels / channel names in the stream (%i).',size(chunk,1),size(stream.buffer,1)); end
        if ~isfield(stream,'marker_pos')
            error('The stream %s has no marker_pos field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        
        % validate marker-related properties
        if exist('markers','var')
            % validate markers argument
            if ~(isequal(markers,[]) || isstruct(markers))
                error('The given Markers argument must be either [] or a struct array, but was: %s',hlp_tostring(markers,10000)); end
            if ~isempty(markers) 
                if ~isfield(markers,'type')
                    error('The given Markers array is lacking the required .type field.'); end
                if ~isfield(markers,'latency')
                    error('The given Markers array is lacking the required .latency field.'); end
            end
            % validate stream fields
            if ~isfield(stream,'mmax')
                error('The stream %s has no mmax field. It might have been overwritten by another function or you use the wrong stream name.',name); end
            if ~isfield(stream,'marker_buffer')
                error('The stream %s has no marker_buffer field. It might have been overwritten by another function or you use the wrong stream name.',name); end
            if ~isfield(stream,'marker_buffer_len')
                error('The stream %s has no marker_buffer_len field. It might have been overwritten by another function or you use the wrong stream name.',name); end            
        end
        
        % validate timestamp-related properties
        if exist('stamp','var')
            % validate stamp argument
            if ~isnumeric(stamp) || ~isscalar(stamp)
                error('The given Timestamp argument must be a scalar numeric value, but was: %s',hlp_tostring(stamp,10000)); end
            % validate stream fields
            if ~isfield(stream,'tmax')
                error('The stream %s has no tmax field. It might have been overwritten by another function or you use the wrong stream name.',name); end
            if ~isfield(stream,'timestamps')
                error('The stream %s has no timestamps field. It might have been overwritten by another function or you use the wrong stream name.',name); end
            if ~isfield(stream,'timestamps_len')
                error('The stream %s has no timestamps_len field. It might have been overwritten by another function or you use the wrong stream name.',name); end
            if ~isfield(stream,'srate')
                error('The stream %s has no srate field. It might have been overwritten by another function or you use the wrong stream name.',name); end
        end
        % unknown error
        rethrow(hlp_handleerror(e));
    catch e
        disp(sprintf('onl_append error: %s',e.message));
    end
end


function perform_update(streamname,dataname,datavalue,code)
persistent cannot_script;
% put the data into the workspace
assignin('base',[streamname '_' dataname],datavalue);
if cannot_script
    % if we previously determined that we cannot work with scripts for some reason,
    % we run the code directly (at some performance cost)
    evalin('base',strrep(sprintf(code),'X',streamname));
else
    try
        % try to run an auto-generated script for this stream and data value (saving interpreter
        % time over running the code directly)
        scriptname = ['update__' streamname '_' dataname];
        evalin('base',scriptname);
    catch e
        code = strrep(sprintf(code),'X',streamname);
        if ~exist(scriptname,'file')
            % try to create the script if it doesn't exist
            try                
                filename = env_translatepath(['functions:/temp/' scriptname '.m']);
                f = fopen(filename,'w+');
                fprintf(f,code);
                fclose(f);
                rehash;
            catch e
                cannot_script = true;
                warning('BCILAB:cannot_create_script','Could not create script "%s" (error message: %s); please make sure that you have write permission in that folder.',filename,e.message);
                evalin('base',code);
                return;
            end
            % then try to run the newly created script
            try
                evalin('base',scriptname);
            catch e
                if strcmp(e.identifier,'MATLAB:UndefinedFunction')
                    cannot_script = true;
                    evalin('base',code);
                else
                    rethrow(e);
                end
            end
        else
            rethrow(e);
        end
    end
end
