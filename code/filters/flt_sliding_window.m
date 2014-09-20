function [signal,state] = flt_sliding_window(varargin)
% Buffer samples and emit overlapping windows.
% function [Signal,State] = flt_sliding_window(Signal,WindowLength,WindowStep)
%
% This filter does *not* generate event-locked segments (see set_makepos for that); instead, it 
% generates *successive segments* with configurable overlap.
%
% This filter turns continuous-time data into epoched data, and handles all time-series fields 
% in the data. It assumes that each time-series field has the same sampling rate.
%
% In:   
%   Signal      : EEGLAB data set, either continuous or epoched
%
%   WindowLength : Window length to emit. (default: 1)
%
%   WindowStep : Step size between successive windows. (default: 0.1)
%
%   TimeUnit : Unit of time parameters. This applies to WindowLength and WindowStep. 
%              Can be 'seconds' or 'samples'. (default: 'seconds')
%
%   State : input state
%
% Out:
%   Signal : segmented signal
%
%   State : output state
%
% Examples:
%   % use default settings
%   eeg = flt_sliding_window(eeg)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-02-12

if ~exp_beginfun('filter') return; end

declare_properties('name','SlidingWindow', 'follows','flt_selchans', 'independent_channels',true, 'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'window_len','WindowLength'}, 1, [], 'Window length to emit.'),...
    arg({'window_step','WindowStep'}, 0.1, [], 'Step between successive windows. If this is set to zero, it is taken to be the same as window_len (i.e., no overlap).'),...
    arg({'time_unit','TimeUnit'}, 'seconds', {'seconds','samples'}, 'Unit of time parameters. This applies to WindowLength and WindowStep.'), ...    
    arg({'handle_events','HandleEvents'}, false, [], 'Handle events. This will segment events in accordance with the chunking of the .data time-series field.'), ...
    arg_sub({'online_options','OnlineOptions'},{},{ ...
        arg({'only_latest_segment','OnlyLatestSegment'},true,[],'Retain only the latest segment. In this case all but the latest segment are dropped from the output during online processing.'), ...    
        arg({'always_at_end','AlwaysAtEnd'},false,[],'Segment always ends on last data sample. In this case the returned segment always ends on the last sample of the data and does not respect WindowStep.'), ...
    },'Online processing options. These only apply when the filter is called from the online system (onl_ functions).'), ...
    arg_nogui({'state','State'}));

if window_step == 0
    window_step = window_len; end

% check for event handling
if handle_events && ~isempty(signal.event)
    disp_once('Warning: the event handling is not yet implemented in flt_segmentation. Dropping events.'); 
    signal.event = [];
end

% perform unit conversion
if strcmp(time_unit,'seconds')
    window_len = max(1,round(window_len*signal.srate));
    window_step = max(1,round(window_step*signal.srate));
elseif ~strcmp(time_unit,'samples')
    error('The time unit must be either samples or seconds, but was: %s',hlp_tostring(time_unit,100));
end

if isempty(state)
    state.buffers = struct();   % data held in buffer
    state.buffer_end_idx = 0;   % cumulative index of the last sample in buffer (1 = first sample in stream)
    state.next_chunk_idx = 1;   % start sample index of the next chunk to be emitted (1 = first sample in stream)
end

signal.pnts = 0;

% are we online and want to return a segment that ends on the last sample?
if online_options.always_at_end && onl_isonline
    
    % for each time-series fields...
    for f = utl_timeseries_fields(signal)
        field = f{1};
        if ~isempty(signal.(field))
            % init state if necessary
            if ~isfield(state.buffers,field)
                state.buffers.(field) = []; end
            
            % concat signal and buffer, cut excess data
            state.buffers.(field) = cat(2,state.buffers.(field),signal.(field));
            state.buffers.(f{1}) = state.buffers.(f{1})(:,max(1,end-window_len+1):end,:,:,:,:,:,:);
            
            if size(state.buffers.(field),2) >= window_len
                % long enough: return data
                signal.(field) = state.buffers.(field);
            else
                % too short: return no data
                signal.(field) = [];
            end
        end
    end
    
else
    % otherwise we return chunks at fixed rate

    % for each time series field...
    for f = utl_timeseries_fields(signal)
        field = f{1};
        if ~isempty(signal.(field))

            % skip already epoched fields
            if size(signal.(field),3) > 1 || ~isempty(signal.epoch)
                disp_once('Note: the time-series field %s is already epoched and will not be further segmented.',field); 
                continue;
            end

            % init state if necessary
            if ~isfield(state.buffers,field)
                state.buffers.(field) = []; end

            % concat data and buffer
            state.buffers.(field) = cat(2,state.buffers.(field),signal.(field));
            state.buffer_end_idx = state.buffer_end_idx + size(signal.(field),2);

            chunks = {};

            % generate chunks and consume buffer
            % while we can emit another chunk...
            while state.next_chunk_idx + window_len-1 <= state.buffer_end_idx
                buffer_start_idx = state.buffer_end_idx-size(state.buffers.(field),2)+1;
                range_in_buffer = (state.next_chunk_idx - buffer_start_idx + 1) + (0:window_len-1);
                chunks{end+1} = state.buffers.(field)(:,range_in_buffer,1,:,:,:,:,:);
                state.next_chunk_idx = state.next_chunk_idx + window_step-1;
            end

            % optionally drop all but latest segment in online mode
            if online_options.only_latest_segment && length(chunks)>1 && onl_isonline
                chunks = chunks(end); end

            % turn into epoched field
            signal.(field) = cat(3,chunks{:});
            signal.pnts = size(signal.(field),2);

            % consume buffer content that we don't need any more
            buffer_start_idx = state.buffer_end_idx-size(state.buffers.(field),2)+1;
            if state.next_chunk_idx > buffer_start_idx
                state.buffers.(field)(:,1:(state.next_chunk_idx - buffer_start_idx),:,:,:,:,:,:) = [];            
                buffer_start_idx = state.next_chunk_idx;
                state.buffer_end_idx = buffer_start_idx + size(state.buffers.(field),2)-1;
            end
        end
    end

end

% update meta-data
if ~isempty(signal.data)
    metadata_field = 'data';
else
    metadata_field = utl_timeseries_fields(signal);
    metadata_field = metadata_field{1};
end
[signal.nbchan,signal.pnts,signal.trials,extra_dims] = size(signal.(metadata_field)); %#ok<NASGU>
signal.xmin = signal.xmax-(signal.pnts-1)/signal.srate;

exp_endfun;

