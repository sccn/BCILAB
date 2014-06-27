function chunk = utl_buffer(chunk,buffer,desired_length)
% Append a chunk to a buffer and trim the buffer to the desired length.
% Buffer = utl_buffer(Chunk,Buffer,DesiredLength)
%
% In:
%   Chunk : a chunk of data (continuous EEGLAB dataset struct)
%
%   Buffer : a buffer of previous data (continuous EEGLAB dataset struct)
%            each non-empty time-series field present in the buffer must have the same number of
%            channels as the corresponding field in the chunk
%
%   DesiredLength : the number of samples of desired output
%
% Out:
%   Buffer : the updated and trimmed buffer
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-12

% input validation
utl_check_fields(chunk,{'data','event','xmax','srate'},'chunk','signal');
utl_check_fields(buffer,{'data','event'},'buffer','signal');
if ~isempty(chunk.event)
    if ~isfield(chunk.event,'latency')
        error('The chunk is missing the required .event.latency field.'); end
    latency_numels = cellfun('prodofsize',{chunk.event.latency});
    if any(latency_numels == 0)
        error('The given chunk has one or more events with empty .latency field. This is not permitted.');
    elseif any(latency_numels ~= 1)
        error('The given chunk has one or more events with a .latency value that is not a scalar. This is not permitted.');
    end
end

% concatenate markers if necessary
if (desired_length ~= size(chunk.data,2)) && (~isempty(buffer.event) || ~isempty(chunk.event))
    if size(chunk.data,2) > desired_length
        if ~isempty(chunk.event)
            % shift the latency of the markers based on how many samples we have to cut
            latency = [chunk.event.latency] - (size(chunk.data,2)-desired_length);
            [chunk.event.latency] = arraydeal(latency);
            % remove excess markers
            chunk.event(latency<0.5) = [];
        end
    else
        if ~isempty(chunk.event)
            % shift the latency of the chunk's markers based on how many samples we have to prepend
            [chunk.event.latency] = arraydeal([chunk.event.latency]+min(size(buffer.data,2),(desired_length-size(chunk.data,2))));
        end
        if ~isempty(buffer.event)
            % shift the latency of the buffer's markers based on how many samples we cut from that
            latency = [buffer.event.latency] - max(0,(size(buffer.data,2)+size(chunk.data,2)-desired_length));
            [buffer.event.latency] = arraydeal(latency);
            % concatenate markers
            if isempty(chunk.event)
                chunk.event = buffer.event(latency>=0.5);
            else
                try
                    chunk.event = [buffer.event(latency>=0.5) chunk.event];
                catch e
                    fprintf('WARNING: Trying to concatenate incompatible event structs in the buffer (%s) and the incoming chunk (%s); dropping buffer contents; error message: %s\n',hlp_tostring(fieldnames(buffer.event)),hlp_tostring(fieldnames(chunk.event)),e.message);
                end
            end
        end
    end
end

% for each time-series field in the chunk...
for fld = utl_timeseries_fields(chunk)
    field = fld{1};    
    
    % make sure that we have the requested field in the previous buffer to avoid special cases below
    if ~isfield(buffer,field)
        buffer.(field) = []; end
    
    % if some other amount of data than what's in the chunk was requested
    if desired_length ~= size(chunk.(field),2)
        if size(chunk.(field),2) < desired_length
            try
                if size(buffer.(field),2) == desired_length
                    % we can do an in-place update without reallocations
                    chunk.(field) = cat(2,buffer.(field)(:,size(chunk.(field),2)+1:end,:,:,:,:,:,:),chunk.(field));
                else
                    % append new samples & cut excess data
                    chunk.(field) = cat(2,buffer.(field),chunk.(field));
                    if size(chunk.(field),2) > desired_length
                        chunk.(field) = chunk.(field)(:,end-desired_length+1:end,:,:,:,:,:,:); end
                end
            catch e
                error('Error trying to concatenate time-series field .%s of buffer (size=%s) and incoming chunk (size=%s); error message: %s.',field,hlp_tostring(size(buffer.(field))),hlp_tostring(size(chunk.field)),e.message);
            end
        else
            % if the chunk is longer than what's requested cut the excess data
            chunk.(field) = chunk.(field)(:,end-desired_length+1:end,:,:,:,:,:,:);
        end
    end
end

% update time-related meta-data
chunk.pnts = size(chunk.data,2);
chunk.xmin = chunk.xmax - (chunk.pnts-1)/chunk.srate;
