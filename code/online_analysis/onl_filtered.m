function [chunk,pipeline] = onl_filtered(pipeline,desired_length,suppress_output,set_online_scope)
% Obtain processed data from a filter pipeline online.
%
% This function returns a chunk of most recent filtered output from a filter pipeline.
% 
% A filter pipeline is a recursive data structure (like a tree) whose nodes are the filter stages and 
% whose edges represent the output of one stage going into the input of another stage. The leaf
% nodes refer to raw online streams (structs in the workspace) and are queried via onl_peek, all other
% nodes are evaluated by calling the filter function on its input data (recursively).
%
% In:
%   Pipeline : previous filter pipeline struct
%
%   DesiredLength : number of samples to get (or 0 to get all new samples)
%
%   SuppressOutput : suppress console output (default: true)
%
%   SetOnlineScope : set the regular online-processing scope (can be turned off if that scope is
%                    already set for some reason) (default: true)
%
% Out:
%   Chunk : EEGLAB dataset struct representing the desired data chunk
%           Can be shorter than desired length if not enough data is available yet
%
%   Pipeline : updated filter pipeline struct
%
%
% Example:
%   % load calibration set
%   raw = io_loadset('calib.set')
%
%   % apply a series of filter to it (the processed set now has a filter expression and initial state)
%   processed = exp_eval(flt_iir(flt_resample(raw,128),[0.5 1],'highpass'));
%
%   % start streaming some data
%   run_readdataset('mystream','action.set');
%   % and put a pipeline on top of it that replicates the processing applied to processed and continues it on new data
%   pip = onl_newpipeline(processed,{'mystream'});
%
%   while 1
%      % generate a 200-sample view into the processed stream
%      [EEG,pip] = onl_filtered(pip,200);
%   end
%
% See also:
%   onl_newpipeline, onl_newstream, onl_append, onl_peek
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-05-13

% check inputs
if nargin < 4
    set_online_scope = true;
    if nargin < 1
        error('You need to pass a pipeline structure.'); end
    if nargin < 2
        error('You need to pass the length of the view that should be generated, in samples.'); end
    if nargin < 3
        suppress_output = true; end
end


% run update_pipeline() with appropriate options
if ~suppress_output
    if ~set_online_scope
        [chunk,pipeline] = update_pipeline(pipeline,desired_length); 
    else
        [chunk,pipeline] = hlp_scope({'disable_expressions',1,'is_online',1},@update_pipeline,pipeline,desired_length); 
    end
else
    if ~set_online_scope
        [console_output,chunk,pipeline] = evalc('update_pipeline(pipeline,desired_length)'); %#ok<ASGLU>
    else
        [console_output,chunk,pipeline] = evalc('hlp_scope({''disable_expressions'',1,''is_online'',1},@update_pipeline,pipeline,desired_length)'); %#ok<ASGLU>
    end
end


function [chunk,p] = update_pipeline(p,desired_length)
% Upate the given filter pipeline and get a chunk most recent output of desired length
% [Chunk,Pipeline] = update_pipeline(Pipeline,DesiredLength)
%
% A pipeline is a recursive data structure (like a tree) whose nodes are the filter stages and 
% whose edges represent the output of one stage going into the input of another stage. The leaf
% nodes refer to raw data streams (structs in the workspace) and are queried via onl_peek, all other
% nodes are evaluated by calling the filter function on its input data (recursively).
%
% At each node we store the filter function (.head) and its arguments (.parts), some of which may be
% dependent filter pipelines themselves, plus some miscellaneous book-keeping data, like a buffer of
% outputs from the last update (.buffer), the previous filter state (.state), and some meta-data
% that characterizes the behavior of the node (.israw, .subrequests, .pipeline_indices).
%
% In:
%   Pipeline : previous filter pipeline struct
%
%   DesiredLength : amount of most recent filtered data to get (or 0 to get all updated data)
%
% Out:
%   Chunk : EEGLAB dataset struct representing the desired data chunk
%           Can be shorter than desired length if not enough data is available yet
%
%   Pipeline : updated filter pipeline struct

% first create potentially missing fields in the pipeline struct
if ~isfield(p,'israw')
    
    % the .israw field encodes whether this stage represents the output of a raw data stream
    p.israw = strcmp(char(p.head),'rawdata');
    
    % the .buffer field contains an EEGLAB-style dataset struct where we buffer the output signal
    p.buffer = struct('data',{[]}, 'smax',{0});

    if ~p.israw
        % the .pipeline_indices field stores which of the input arguments are pipelines (that we need to
        % update recursively)
        p.pipeline_indices = find(cellfun(@(x)all(isfield(x,{'head','parts'})),p.parts));
        
        % the .stateful field says whether the filter function returns state
        if ~isfield(p,'stateful')
            p.stateful = p.israw || is_stateful(p.head); end
        
        % the .state field contains the previous state of the filter function
        if p.stateful && ~isfield(p,'state')
            p.state = []; end
        
        % the .subrequests field describes how much data is requested by this stage from each of its
        % input pipeline stages (this depends on DesiredLength, which we needs to remain constant)
        if ~isfield(p,'subrequests')
            p.subrequests = nan(1,length(p.pipeline_indices));
        elseif length(p.subrequests) ~= length(p.pipeline_indices)
            warn_once('BCILAB:onl_predict:inconsistent_pipeline','A filter pipeline node was encountered with a .subrequests field that does not match its sub-pipelines.');
            p.subrequests = nan(1,length(p.pipeline_indices));
        end
        
        % fill in any undetermined request lengths: if we are stateful, we will request 0, which means
        % "all new data", and if we are stateless, we request the amount of data that is being requested
        % from us (assuming that this pipeline stage is not a rate-changing stateless (i.e. epoch-based)
        % filter). For any such filter to work, subrequests must have already been initialized properly
        % in utl_add_online.
        p.subrequests(isnan(p.subrequests)) = desired_length * ~p.stateful;
    end
end

% get the desired data
if ~p.israw
    % reading from a regular filter stage

    % update the inputs to the current stage (recursively)
    inputs = p.parts;
    for n = 1:length(p.pipeline_indices)
        k = p.pipeline_indices(n);
        % update input pipelines and replace the corresponding input signals by their buffers
        [inputs{k},p.parts{k}] = update_pipeline(p.parts{k},p.subrequests(n));
    end

    % get a chunk of output by applying the filter for this pipeline stage to the inputs
    if p.stateful
        % call stateful filter function (appending previous state, retrieving new state)
        [chunk,p.state] = p.head(inputs{:}, 'state',p.state);
        % increment .smax with what we got (p.head might change the sampling rate)
        chunk.smax = p.buffer.smax + size(chunk.data,2);
    else
        % call stateless filter function on the input arguments
        chunk = p.head(inputs{:});
        % the closest to an .smax in a stateless filter is the .smax of its inputs (we assume that it
        % does not change the rate)
        chunk.smax = inputs{p.pipeline_indices(1)}.smax;
    end
    
    if desired_length ~= 0
        % if any other amount of data than the most recent chunk was requested
        % as output, concatenate with previous p.buffer contents and trim to desired size

        % for each time-series field in the chunk...
        for fld = utl_timeseries_fields(chunk)
            field = fld{1};
            
            % make sure that we have the requested field in the previous buffer to avoid special cases below
            if ~isfield(p.buffer,field)
                p.buffer.(field) = []; end

            % if some other amount of data than what's in the chunk was requested
            if desired_length ~= size(chunk.(field),2)
                if size(chunk.(field),2) < desired_length
                    % stateful filters must concatenate the newly-produced outputs with the
                    % data that they already have in the buffer from previous invocations
                    if p.israw || p.stateful
                        if size(p.buffer.(field),2) == desired_length
                            % we can do an in-place update without reallocations
                            % Note: if you get an error here, a filter (p.head) likely had returned different
                            %       numbers of channels across invocations.
                            chunk.(field) = [p.buffer.(field)(:,size(chunk.(field),2)+1:end) chunk.(field)];
                        else
                            % append new samples & cut excess data
                            % Note: if you get an error here, a filter (p.head) likely had returned different
                            %       numbers of channels across invocations.
                            chunk.(field) = [p.buffer.(field) chunk.(field)];
                            if size(chunk.(field),2) > desired_length
                                chunk.(field) = chunk.(field)(:,end-desired_length+1:end); end
                        end
                    end
                else
                    % if the chunk is longer than what's requested cut the excess data
                    chunk.(field) = chunk.(field)(:,end-desired_length+1:end);
                end
            end        
        end
    end
    
    % ensure consistency of time-related meta-data
    chunk.pnts = size(chunk.data,2);
    try
        chunk.xmin = chunk.xmax - (chunk.pnts-1)/chunk.srate;
    catch
        if ~isfield(chunk,'xmin')
            chunk.xmin = 0; end
        chunk.xmax = chunk.xmin + (chunk.pnts-1)/chunk.srate;
    end
else
    % read from a raw stream: p.parts of the stage holds {stream_name,channel_range)
    if desired_length ~= 0
        % get the requested number of samples
        chunk = onl_peek(p.parts{1},desired_length,'samples',p.parts{2});
    else
        % get the most recent samples since our buffer's smax
        chunk = onl_peek(p.parts{1},p.buffer.smax,'index',p.parts{2});
    end
end

% replace old buffer
p.buffer = chunk;
