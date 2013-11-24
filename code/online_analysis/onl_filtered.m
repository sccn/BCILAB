function [chunk,p] = onl_filtered(p,desired_length,suppress_output,set_online_scope)
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
%   DesiredLength : number of samples to get (or 0 to get all new samples) (default: 0)
%
%   SuppressOutput : suppress console output (default: true)
%
%   SetOnlineScope : set the regular online-processing scope (can be turned off for efficiency if
%                    that scope is already set for some reason) (default: true)
%
% Out:
%   Chunk : EEGLAB dataset struct representing the desired data chunk
%           Can be shorter than desired length if not enough data is available at the moment; if the
%           chunk is epoched, the desired length is ignored.
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
        desired_length = 0; end
    if nargin < 3
        suppress_output = true; end
end

% run update_pipeline() with appropriate options
if suppress_output
    if set_online_scope
        [console_output,chunk,p] = evalc('hlp_scope({''disable_expressions'',1,''is_online'',1},@update_pipeline,pipeline)'); %#ok<ASGLU>
    else
        [console_output,chunk,p] = evalc('update_pipeline(pipeline)'); %#ok<ASGLU>
    end
else
    if set_online_scope
        [chunk,p] = hlp_scope({'disable_expressions',1,'is_online',1},@update_pipeline,p); 
    else
        [chunk,p] = update_pipeline(p); 
    end
end

% if a desired length was specified and the chunk is not epoched
if desired_length && ~isstruct(chunk.epoch)    
    try
        % concatenate the chunk with previous buffer contents and trim to desired size
        p.buffer = utl_buffer(chunk,p.buffer,desired_length);
    catch e
        % handle missing .buffer field (first-time use)
        if ~isfield(p,'buffer')        
            p.buffer = utl_buffer(chunk,struct('data',[],'event',[]),desired_length);
        else
            rethow(e);
        end
    end
    chunk = p.buffer;
end


function [chunk,p] = update_pipeline(p)
% Update the given filter pipeline and get a chunk of the newly appended output
% [Chunk,Pipeline] = update_pipeline(Pipeline)
%
% A pipeline is a recursive data structure (like a tree) whose nodes are the filter stages and 
% whose edges represent the output of one stage going into the input of another stage. The leaf
% nodes refer to raw data streams (structs in the workspace) and are queried via onl_peek, all other
% nodes are evaluated by calling the filter function on its input data (recursively).
%
% At each node we store the filter function (.head) and its arguments (.parts), some of which may be
% input filter pipelines themselves, plus some miscellaneous book-keeping data. These include:
%  * .israw : true for nodes that represent raw data
%  * .pipeline_indices : indices of those input arguments that are pipelines themselves (if not raw)
%  * .stateful : true if the node has state
%  * .state : previous filter state, if stateful
%  * .smax : number of samples seen so far (if israw)
%
% In:
%   Pipeline : previous filter pipeline struct
%
% Out:
%   Chunk : EEGLAB dataset struct representing newly appended data, filtered
%
%   Pipeline : updated filter pipeline struct

inputs = p.parts;
if p.subnodes
    % update input pipelines to the current node and store the results in inputs
    for k=p.subnodes
        [inputs{k},p.parts{k}] = update_pipeline(p.parts{k}); end
    % process the inputs by calling the respective filter function
    if p.stateful
        [chunk,p.state] = p.head(inputs{:},'state',p.state);
    else
        chunk = p.head(inputs{:});
    end
else
    % get the most recent samples since our buffer's smax from a raw stream: 
    % inputs holds the cell array {stream_name,channel_range)
    chunk = onl_peek(inputs{1},p.smax,'index',inputs{2});
    p.smax = chunk.smax;
end
