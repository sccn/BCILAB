function onl_append(name, chunk, stamp)
% Append a block of raw data to a stream.
% onl_append(Name, Chunk, Timestamp)
%
% The number of columns in the appended data block must match the number of channels in the stream,
% and the data must have been sampled at the rate specified in the stream's .srate field.
% Predictions are always being made with respect to the most recently appended sample.
%
% The sequence of samples fed into a stream (over the course of multiple calls to onl_append) must
% have no omissions (i.e., samples must not be skipped) or repetitions (i.e., samples that have been
% supplied in a previous call must not be fed again). If samples are not available from the data
% source or were missed (for whatever reason), the corresponding number of zero-valued samples
% should be fed).
%
% Processing (i.e., appending samples and calculating predictions from them) does not have to happen
% in real time (i.e., the actual wall-clock time does not affect the calulated results).
%
% In:
%   Name : variable name of the stream to which the data should be appended
%
%   Chunk : [#Channels x #Samples] matrix of raw data
%
%   Timestamp : optional time stamp for the last sample in the chunk, in seconds, in some arbitary
%               time domain which is consistent across all streams used by dependent predictors
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
% See also:
%   onl_newstream, onl_newpredictor, onl_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03


try
    if nargin <= 2 || isempty(stamp)
        
        try
            
            % put the chunk into the workspace
            assignin('base',[name '_chunk'],chunk);            
            try
                % append the chunk to the stream (atomically); we do this by running an 
                % auto-generated mini script, to save the interpreter time
                scriptname = ['update__' name];
                evalin('base',scriptname);
            catch e
                if isempty(chunk)
                    return; end;
                
                % check what's wrong...
                persistent cannot_script; %#ok<TLEV>
                if cannot_script
                    % if we previously determined that we cannot work with scripts for some reason, 
                    % we run the code directly (at some performance cost)
                    evalin('base',['[' name '.smax,' name '.buffer(:,1+mod(' name '.smax:' name '.smax+size(' name '_chunk,2)-1,' name '.buffer_len))] = deal(' name '.smax + size(' name '_chunk,2),' name '_chunk);']);
                else
                    % does the script exist?
                    if ~exist(scriptname,'file')
                        try
                            % try to create it (likely a first-time use)
                            filename = env_translatepath(['functions:/temp/' scriptname '.m']);
                            f = fopen(filename,'w+');
                            fprintf(f,['[' name '.smax,' name '.buffer(:,1+mod(' name '.smax:' name '.smax+size(' name '_chunk,2)-1,' name '.buffer_len))] = deal(' name '.smax + size(' name '_chunk,2),' name '_chunk);']);
                            fclose(f);
                            rehash;
                        catch
                            cannot_script = true;
                            % error during creation
                            warning('BCILAB:cannot_create_script','Could not create script "%s"; please make sure that you have write permission in that folder.',filename);
                            error('scripting problem...');
                        end
                        % now re-run it
                        evalin('base',scriptname);
                    else
                        % script did exist: error likely had a different cause (caught below...)
                        rethrow(e);
                    end
                end
            end
        catch le
            
            % nothing to do?
            if isempty(chunk)
                return; end;
            
            % error: interpret the source of the error & display an appropriate message
            try
                stream = evalin('base',name);
            catch
                error(['A stream named ' name ' does not exist.']);
            end
            if ndims(chunk) ~= 2
                error('onl_append expects a [#Channels x #Samples] matrix.');
            elseif ~isstruct(stream) || ~isscalar(stream)
                error(['The stream ' name ' is not a 1x1 struct. Probably it was accidentally overwritten by another function.']);
            elseif ~isfield(stream,'buffer')
                error(['The stream ' name ' has no buffer. Probably it was accidentally overwritten by another function.']);
            elseif ~isnumeric(stream.buffer)
                error(['The stream ' name ' has an invalid buffer. Probably it was accidentally overwritten by another function.']);
            elseif ~isfield(stream,'buffer_len')
                error(['The stream ' name ' has no buffer_len field. Probably it was accidentally overwritten by another function.']);
            elseif ~isscalar(stream.buffer_len) || ~isnumeric(stream.buffer_len)
                error(['The stream ' name ' has an invalid buffer_len field. Probably it was accidentally overwritten by another function.']);
            elseif size(stream.buffer,2) ~= stream.buffer_len
                error(['The stream ' name ' has a buffer with an invalid length. Probably it was accidentally overwritten by another function.']);
            elseif ~isfield(stream,'smax')
                error(['The stream ' name ' has no smax field. Probably it was accidentally overwritten by another function.']);
            elseif ~isscalar(stream.smax) || ~isnumeric(stream.smax)
                error(['The stream ' name ' has an invalid smax field. Probably it was accidentally overwritten by another function.']);
            elseif size(stream.buffer,1) ~= size(chunk,1)
                error('Number of channels in the supplied chunk does not match the number of channels / channel names in the stream.');
            elseif strcmp('MATLAB:UndefinedFunction',le.identifier)
                cannot_script = true;
                error('BCILAB:script_failure','Execution of script "%s" failed; falling back to interpreted code',scriptname);
            else
                rethrow(le);
            end
            
        end
        
    else
        % time stamp was supplied, need to do non-trivial time-stamp updating: get the stream
        stream = evalin('base',name);
        
        % update the buffer
        stream.buffer(:,1+mod(stream.smax:stream.smax+size(lastchunk,2)-1,stream.buffer_len)) = chunk;
        stream.smax = stream.smax + size(lastchunk,2);
        % insert timestamp at the end if the buffer (contains the time in timestamp domain and in
        % 0-based stream time domain)
        stream.timestamp_ptr = 1+mod(stream.timestamp_ptr,stream.timestamp_avgs);
        stream.timestamps(stream.timestamp_ptr,:) = [stamp stream.smax/stream.srate];
        % estimate time lag (offset of time stamp time relative to 0-based stream time)
        lag = sum(stream.timestamps(:,1) - stream.timestamps(:,2))/nnz(stream.timestamps(:,2));
        % derive xmin/xmax in the time stamp domain
        stream.xmax = stream.smax/stream.srate + lag;
        stream.xmin = (1+max(0,stream.smax-stream.buffer_len))/stream.srate + lag;
        
        % atomically write back the updated stream
        assignin('base',name,stream);
    end
catch e
    disp('An exception occured during onl_append:');
    env_handleerror(e);
end
