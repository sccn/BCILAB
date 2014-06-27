function pip = utl_resolve_streams(pip,streams,chanlocs)
% resolve the stream that each rawdata node in the given filter chain requires, given a list
% of workspace stream names (or stream structs), and a subset of channels needed by above pipelines (where empty means 'all')
%
% In:
%   Pipeline : a filter chain or filter graph (cell array of filter chains) with rawdata nodes to resolve
%
%   Streams : a cell array of workspace stream names, or a cell array of stream/dataset structs,
%             or a stream bundle (struct with .streams{} cell array); these are the candidates
%             that shall be subsituted into the open rawdata nodes
%
%   Chanlocs : optionally the chanlocs structure that is expected at the outlet of the filter
%              chain, as cell array of chanlocs structs in case the pipeline is a filter graph; if 
%              known, it can help in situations where a candidate stream provides fewer channels than
%              technically required by the filter chain -- but not effectively required as the chain
%              prunes some of them further down the processing
%
% Out:
%   Pipeline : the pipeline with matched streams filled in at the @rawdata nodes;
%              if streams was a name array, then the rawdata nodes are augmented with the stream name
%              (for online processing); otherwise, the expression of the respective stream is
%              substituted in place of the respective rawdata node
%
% See also:
%   bci_train, bci_predict, onl_newpredictor
%
% Notes:
%   The matching is generally case insensitive (both for channel labels and channel types).
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-23

% turn streams into a cell array
if isstruct(streams) && isscalar(streams) && isfield(streams,'streams')    
    streams = streams.streams; end
if ~iscell(streams)
    streams = {streams}; end

% resolve stream properties
for s = 1:length(streams)
    if ischar(streams{s})
        % stream given as a workspace variable name:
        % - can read out the chanlocs.labels directly from it
        % - we take the types as declared by the stream
        % - we substitute the variable name into the rawdata node, as we will refer to this stream online
        tmp = evalin('base',streams{s});
        streams{s} = struct('labels',{{tmp.chanlocs.labels}},'types',{tmp.types},'substitute',streams{s});
    elseif all(isfield(streams{s},{'chanlocs','tracking'}))
        % stream is given as a full data set:
        % - can read out the chanlocs & their types directly from it
        % - we substitute the streams expression in there
        streams{s} = struct('labels',{{streams{s}.chanlocs.labels}},'types',{unique({streams{s}.chanlocs.type})},'substitute',streams{s}.tracking.expression);
    elseif all(isfield(streams{s},{'head','parts'}))
        % stream is an unevaluated expression:
        % - we try to infer the chanlocs & their types from it without evaluating too much of it
        % - we substitute the stream expression itself
        streams{s} = resolve_stream_properties(streams{s});
    end
end

% now resolve rawdata nodes using the stream information
if ~iscell(pip)
    % pipeline given as single filter chain
    if exist('chanlocs','var') && ~isempty(chanlocs)
        labels = {chanlocs.labels};
        if isfield(chanlocs,'type')
            types = unique({chanlocs.type});
        else
            types = [];
        end
    else
        labels = [];
        types = [];
    end
    pip = resolve_rawdata(pip,streams,labels,types);
else
    % pipeline given as a cell array: resolve each one by one
    if exist('chanlocs','var') && ~isempty(chanlocs)
        if ~iscell(chanlocs)
            error('If given, chanlocs should be a cell array with one cell for each chain in the filter graph.'); end
        if length(chanlocs) ~= length(pip)
            error('The Chanlocs array has a different number of elements than the filter graph.'); end
        for p = 1:length(pip)
            labels{p} = {chanlocs{p}.labels}; %#ok<AGROW>
            if isfield(chanlocs{p},'type')
                types{p} = unique({chanlocs{p}.type}); %#ok<AGROW>
            else
                types{p} = [];
            end
        end
    else
        labels = repmat({[]},1,length(pip));
        types = repmat({[]},1,length(pip));
    end
    for p = 1:length(pip)
        pip{p} = resolve_rawdata(pip{p},streams,labels{p},types{p}); end
end


function pip = resolve_rawdata(pip,streams,chan_subset,type_subset)
% Resolve rawdata nodes in a single filter chain into streams.
% 
% In:
%   Pipeline : the filter chain
%
%   Streams : cell array of structs with fields .labels, .types, and .substitute
%             (.labels is a cell array of labels provided by the stream, .types is a cell array of 
%              channel types (unique) provided by the stream, and .substitute is  the object that 
%              should be used for substitution in the respective rawdata slot if it matches)
%
%   RequiredChannels : full set of channel labels that are required from this pipeline (may be a
%                      subset of the channel labels expected at the @rawdata node, in which case 
%                      fewer channels need to be provided by a matching stream); may be a cell array
%                      of channel labels or index range, or [] (which stands for "all channels" of this stage)
%
%   RequiredTypes : full set of channel types that are required from this pipeline (may be a subset
%                   of the types expected at the rawdata node, in which case fewer channeels need to 
%                   be provided by the matching stream); may also be [] (which stands for "all types" of this stage)
%
% Out:
%   Pipeline : the filter chain with @rawdata nodes substituted appropriately
%

phead = char(pip.head);
if ~strcmp(phead,'rawdata')
    
    % --- intermediate node ---

    % figure out the needed channels & types along the way
    if strcmp(phead,'flt_selchans')
        % a selchans node: informs us that we need only a subset of channels from our input streams
        
        % find the last index in the pipeline's arguments that is 'channels'...
        % the channel selection is the argument following that (it is a name-value pair)
        ci = find(strcmp(pip.parts,'channels'),1,'last');
        selected_chans = pip.parts{ci+1};
        % our current channel working subset is either a subset of the selection done by this
        % flt_selchans or refers to "all" channels
        if isempty(chan_subset)
            % "all" channels required: we need the entire chan selection at this level
            chan_subset = selected_chans;
        elseif isnumeric(chan_subset)
            % we need an index subrange of pip.parts
            chan_subset = selected_chans(chan_subset);
        else
            % we known which named subset we need; the flt_selchans does not provide additional info
        end
    elseif strcmp(phead,'flt_seltypes')
        % a seltypes node: informs us that we need only a subset of channel types from here on
        if isempty(type_subset)
            % find the last index in the pipeline's arguments that is 'types'...
            % the channel selection is the argument following that (it is a name-value pair)
            ti = find(strcmp(pip.parts,'chantypes'),1,'last');
            type_subset = pip.parts{ti+1};
            if ~iscell(type_subset)
                type_subset = {type_subset}; end
        else
            % note: if we already have a type subset from farther up in the processing chain we
            %       don't need additional info from this stage
        end
        
    else
        % a generic node: check if this one mixes channels -- in which case we need all channels & types from below
        props = hlp_microcache('props',@arg_report,'properties',pip.head,{hlp_fileinfo(pip.head,[],'hash')});
        if ~props.independent_channels
            % the node doesn't have independent channels: need all channnels & types from below
            chan_subset = [];
            type_subset = [];
        end
    end

    % and follow down its inlets
    for p = find(cellfun(@(p)all(isfield(p,{'head','parts'})),pip.parts))
        pip.parts{p} = resolve_rawdata(pip.parts{p},streams,chan_subset,type_subset); end
    
else
    % --- reached a rawdata node: resolve it using streams ---
    
    % find out what types we need
    if isempty(type_subset)
        % all types from here are needed
        needed_types = pip.parts{2};
    elseif iscell(type_subset)
        % can restrict the types according to the already known type subset (from further up the pipeline)
        needed_types = intersect(pip.parts{2},type_subset);
    else
        error('Unsupported type subset format.');
    end
    
    % find out what channels we need
    if isempty(chan_subset)
        % all channels are needed
        needed_channels = pip.parts{1};
    elseif iscell(chan_subset)
        % a particular named set of channels is needed
        needed_channels = chan_subset;
    elseif isnumeric(chan_subset)
        % an index range referring to the base channel set is needed
        needed_channels = pip.parts{1}(chan_subset);
    else
        error('Unsupported channel subset format.');
    end
    
    satisfies_channels = logical([]); % the stream can provide all the channels required by the pipeline (ideal case if non-ambiguous)
    satisfies_types = logical([]);    % the stream has the same number of channels and set of types as required by the pipeline (happens when the online plugin doesn't supply the correct channel labels...)
    % for each stream...
    for s = 1:length(streams)
        % find out if the stream matches
        satisfies_channels(s) = isempty(fast_setdiff(lower(needed_channels),lower(streams{s}.labels)));
        satisfies_types(s) = length(pip.parts{1}) == length(streams{s}.labels) && isequal(sort(lower(needed_types)),lower(streams{s}.types));
    end
    
    % resolve the stream used for this rawdata node
    if nnz(satisfies_channels) == 1
        % found exactly one applicable stream
        strm = streams{satisfies_channels};
        if ischar(strm.substitute)
            % substitute a variable name into the rawdata node
            pip = struct('head',@rawdata,'parts',{{strm.substitute,set_chanid(strm.labels,needed_channels)}});
        else
            % substitute the rawdata node by the stream's expression
            pip = strm.substitute;
        end
    else
        
        % no clear-cut case: check whether the channel numbers & set of types matches
        if nnz(satisfies_channels) > 1
            % found more than one stream: include the type constraint
            if nnz(satisfies_channels & satisfies_types) == 1
                % now it is non-ambiguous, but warn
                warning('BCILAB:onl_newpredictor:ambiguous','More than one of the supplied/available streams could supply the required channels; falling back to basic channel type requirements to resolve the ambiguity.');
                strm = streams{satisfies_channels & satisfies_types};
                if ischar(strm.substitute)
                    % stream given as a variable name: augment the rawdata node
                    pip = struct('head',@rawdata,'parts',{{strm.substitute,set_chanid(strm.labels,needed_channels)}});
                else
                    % substitute the rawdata node by the stream's expression
                    pip = strm.substitute;
                end
            elseif nnz(satisfies_channels & satisfies_types) > 1
                % still ambiguous
                error('BCILAB:onl_newpredictor:ambiguous','More than one of the supplied/available streams could supply the required channels & types; please provide the desired subset of streams.');
            else
                % the additional constraint eliminates all options
                error('BCILAB:onl_newpredictor:ambiguous','More than one of the supplied/available streams could supply the required channels; please provide the desired subset of streams.');
            end
        elseif nnz(satisfies_types) == 1
            % found exactly one stream which has the correct number of channels and the same set of types
            warning('BCILAB:onl_newpredictor:ambiguous','No stream was found which contains the required channels; however, a stream was found to have the same number of channels and set of types as required by the pipeline.');
            strm = streams{satisfies_types};
            if ischar(strm)
                % stream given as a variable name: augment the rawdata node
                pip = struct('head',@rawdata,'parts',{{strm,1:length(strm.labels)}});
            else
                % substitute the rawdata node by the stream's expression
                pip = strm.substitute;
            end
        else
            % found no applicable stream
            error('BCILAB:onl_newpredictor:nostream',['None of the supplied/available streams has the required channels: ' hlp_tostring(needed_channels) '.']);
        end
        
    end
end

function output = resolve_stream_properties(stream)
% recursively infer the channel labels & types, as well as the substitution expression produced by this pipeline stage
if isfield(stream,'chanlocs')
    % we have reached a data set node: derive all properties
    output.chanlocs = stream.chanlocs;
    output.labels = {stream.chanlocs.labels};
    output.types = unique({stream.chanlocs.type});
    try
        output.substitute = stream.tracking.expression;
    catch %#ok<CTCH>
        error('Reached a dataset node which is lacking an associated expression. This should have been resolved previously using utl_check_dataset or utl_check_bundle.');
    end
elseif all(isfield(stream,{'head','parts'}))
    % we are looking at a symbolic expression
    output.substitute = stream;
    % check if it is a special known one
    shead = char(stream.head);
    switch(shead)
        case 'io_loadset'
            % the chanlocs info can be obtained from the underlying raw data set
            tmp = exp_eval_optimized(stream);
            output.chanlocs = tmp.chanlocs;
            output.labels = {tmp.chanlocs.labels};
            output.types = unique({tmp.chanlocs.type});
        case 'flt_selchans'
            % a channel subset has been selected here: get the chanlocs from further below...
            subprops = resolve_stream_properties(stream.parts{find(cellfun(@(p)all(isfield(p,{'head','parts'})),stream.parts),1)});
            % and get the selected channel labels
            selection = stream.parts{find(strcmp(stream.parts,'channels'),1,'last')+1};
            output.chanlocs = subprops.chanlocs(set_chanid(subprops.chanlocs,selection));
            output.labels = {output.chanlocs.labels};
            output.types = unique({output.chanlocs.type});
        case 'flt_seltypes'
            % a type subset has been selected here: get the chanlocs from further below
            subprops = resolve_stream_properties(stream.parts{find(cellfun(@(p)all(isfield(p,{'head','parts'})),stream.parts),1)});
            % and get the selected type labels
            selection = stream.parts{find(strcmp(stream.parts,'types'),1,'last')+1};
            if ~iscell(selection)
                selection = {selection}; end
            % now restrict both the types and the chanlocs appropriately
            [dummy,ia] = intersect(lower({subprops.chanlocs.type}),lower(selection)); %#ok<ASGLU>
            output.chanlocs = subprops.chanlocs(ia);
            output.labels = {output.chanlocs.labels};
            output.types = unique({output.chanlocs.type});
        case 'exp_block'
            output = resolve_stream_properties(stream.parts{2});
        otherwise
            % generic expression node: check if this node mixes up channels ...
            props = hlp_microcache('props',@arg_report,'properties',stream.head,{hlp_fileinfo(stream.head,[],'hash')});
            if isfield(props,'independent_channels') && props.independent_channels
                % has independent channels: just take the output from below
                output = resolve_stream_properties(stream.parts{find(cellfun(@(p) any(isfield(p,{'parts','chanlocs'})),stream.parts),1)});
                output.substitute = stream;
            else
                % does not have independent channels: need to evaluate the sub-expression
                tmp = exp_eval_optimized(stream);
                output.chanlocs = tmp.chanlocs;
                output.labels = {tmp.chanlocs.labels};
                output.types = unique({tmp.chanlocs.type});
            end
    end
end
