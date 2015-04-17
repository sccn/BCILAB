function pip = utl_resolve_streams(pip,streams,chanlocs)
% Substitute "rawdata" placeholder nodes in a given filter pipeline by given candidate streams.
%
% This function transforms a filter pipeline expression into one where any placeholder nodes that
% represent a raw data stream ("rawdata" nodes) have been resolved into a matching stream data
% structure from a provided set of candidates (Streams argument).
%
% In:
%   Pipeline : a filter chain or filter graph (cell array of filter chains) with "rawdata"
%              placeholder leaf nodes that shall be substituted by given matching candidate streams
%              (note: e.g., the .tracking.online_expression field of a filtered data set is a valid 
%              filter chain)
%
%   Streams : candiate streams to substitute: this is either a cell array of workspace stream names, 
%             or a cell array of stream/dataset structs, or a stream bundle (struct with .streams{} 
%             cell array); these are the candidates that shall be substituted into the rawdata nodes
%
%   Chanlocs : optionally the chanlocs structure that is expected at the outlet of the filter
%              chain, or a cell array of chanlocs structs in case the pipeline is a filter graph; if 
%              known, it can serve to disambiguate situations where a candidate stream provides 
%              fewer channels than technically required by the filter chain while perhaps not all are
%              effectively required due to the processing chain pruning some of them further down
%              in the pipeline
%
% Out:
%   Pipeline : the pipeline with matched streams filled in at the @rawdata nodes;
%              if streams was a name array, then the rawdata nodes are annotated with the stream name
%              (for online processing using onl_*); otherwise, the expression of the respective stream 
%              is substituted in place of the respective rawdata node
%
% See also:
%   bci_train, bci_predict, onl_newpredictor
%
% Notes:
%   The matching of names is generally case insensitive (both for channel labels and channel types).
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-23
dp;

% turn pipeline into a cell array & sanity-check
was_cell = iscell(pip);
if ~was_cell
    pip = {pip}; end
for p=1:length(pip)
    if ~all(isfield(pip{p},{'head','parts'}))
        error('The given filter pipeline must be an expression (struct with fields .head and .parts) or a cell array of such expressions, but was: %s',hlp_tostring(pip,1000)); end
end

% turn streams into a cell array & also handle the special case where a stream bundle is given
if isstruct(streams) && isscalar(streams) && isfield(streams,'streams')    
    streams = streams.streams; end
if ~iscell(streams)
    streams = {streams}; end

% turn chanlocs into a cell array
if nargin < 3
    chanlocs = repmat({[]},1,length(pip)); end
if ~iscell(chanlocs) || iscellstr(chanlocs)
    chanlocs = {chanlocs}; end
if length(chanlocs) ~= length(pip)
    error('If the given pipeline has more than one filter chain (i.e., is a multi-element cell array), then the Chanlocs, if provided, must be a cell array of chanlocs structs with one cell per filter chain.'); end

% determine properties of each provided stream
% (available channel labels, channel types and the data structure to substitute in place of the stream)
for s = 1:length(streams)
    if ischar(streams{s})
        % stream given as a workspace variable name:
        % - can read out the chanlocs.labels directly from it
        % - we take the types as declared by the stream
        % - if matching we annotate the rawdata node with the workspace name, as we will refer to this stream online
        if ~isvarname(streams{s})
            error('If candidate streams are given as strings then each of them must refer to a data structure in the workspace, and therefore be a valid variable name, but was instead: %s',hlp_tostring(streams{s})); end            
        try
            tmp = evalin('base',streams{s});
        catch e
            error('A stream with name "%s" was not found in the base workspace.',streams{s});
        end
        if ~all(isfield(tmp,{'chanlocs','types'}))
            error('The data structure named "%s" in the base workspace is not a valid stream (needs to have both .chanlocs and .types fields), but was: %s.',streams{s},hlp_tostring(streams{s},1000)); end
        streams{s} = struct('labels',{{tmp.chanlocs.labels}},'types',{unique(tmp.types)},'substitute',streams{s});
    else
        % stream is a dataset or dataset expression:
        % - can read out the channel labels and types directly from it
        % - if matching the rawdata node will be substituted by the expression associated with the dataset
        if ~isfield(streams{s},'chanlocs')
            streams{s} = exp_eval_optimized(streams{s}); end
        if all(isfield(streams{s},{'chanlocs','tracking'})) && isfield(streams{s}.tracking,'expression')
            streams{s} = struct('labels',{{streams{s}.chanlocs.labels}},'types',{unique({streams{s}.chanlocs.type})},'substitute',streams{s}.tracking.expression);
        else
            error('The given stream object %s was not recognized (must be either a string or a BCILAB expression that either is a dataset struct or evaluates into one), but was: %s',streams{s},hlp_tostring(streams{s},1000)); 
        end
    end
end

% for each filter chain in the given pipeline...
for p=1:length(pip)
    % get supplementary information from the chanlocs, if any...
    if ~isempty(chanlocs{p}) 
        if iscellstr(chanlocs{p})
            labels = chanlocs{p};
        elseif isfield(chanlocs{p},'labels')
            labels = {chanlocs{p}.labels};
        else
            error('The chanlocs set #%i is neither a cell array of strings nor an eeglab chanlocs struct, but was: %s',p,hlp_tostring(chanlocs{p},1000));
        end
        if isfield(chanlocs{p},'types')
            types = unique({chanlocs.type});
        else
            types = {};
        end
    else
        labels = {};
        types = {};
    end
    % ... and resolve rawdata nodes in the filter chain given candidate streams and supp info
    pip{p} = resolve_rawdata(pip{p},streams,labels,types);
end

if ~was_cell
    pip = pip{1}; end

function pip = resolve_rawdata(pip,candidate_substitutions,required_channels,required_types)
% Resolve rawdata nodes in a single filter chain into streams.
% 
% In:
%   Pipeline : the filter chain
%
%   CandidateSubstitutions : cell array of candidate substitutions; structs with fields .labels,
%              .types, and .substitute (.labels is a cell array of labels provided by the stream,
%              .types is a cell array of channel types (unique) provided by the stream, and
%              .substitute is the object that should be used for substitution in the respective
%              rawdata slot if it matches)
%
%   RequiredChannels : full set of channel labels that are required from this pipeline (may be a
%                      subset of the channel labels expected at the @rawdata node, in which case
%                      fewer channels need to be provided by a matching stream); may be a cell array
%                      of channel labels or index range, or {} (which stands for "all channels" of
%                      this stage)
%
%   RequiredTypes : full set of channel types that are required from this pipeline (may be a subset
%                   of the types expected at the rawdata node, in which case fewer channels need to
%                   be provided by the matching stream); may also be {} (which stands for "all
%                   types" of this stage)
%
% Out:
%   Pipeline : the filter chain with @rawdata nodes substituted appropriately
%

if ~iscellstr(required_channels)
    error('The given set of required channels must be a cell array of strings (possibly empty), but was: %s',hlp_tostring(required_channels,1000)); end
if ~iscellstr(required_types)
    error('The given set of required types must be a cell array of strings (possibly empty), but was: %s',hlp_tostring(required_types,1000)); end

% check the type of node to see if we hit a lead node or whether we can refine the required types
switch char(pip.head)    
    case 'rawdata'
        % reached a rawdata leaf node: determine the best substitution

        % finalize the set of required channel labels needed by the pipeline
        if isempty(required_channels)
            % all channels are needed
            required_channels = pip.parts{1};
        elseif isnumeric(required_channels)
            % an index range referring to the base channel set is needed
            required_channels = pip.parts{1}(required_channels);
        end

        % finalize the set of required channel types needed by the pipeline
        if isempty(required_types)
            % all types from here are needed
            required_types = pip.parts{2};
        else
            % can restrict the types according to the already known type subset (from further up the pipeline)
            required_types = intersect(pip.parts{2},required_types);
        end

        % determine the degree of matching for each stream candidate
        for s = length(candidate_substitutions):-1:1
            % exact matches of the stream's provided channels and the pipeline's input channels
            matches_labels_and_order(s) = isequal(lower(pip.parts{1}),lower(candidate_substitutions{s}.labels));
            matches_labels(s) = isequal(sort(lower(pip.parts{1})),sort(lower(candidate_substitutions{s}.labels)));
            matches_chancount(s) = length(pip.parts{1}) == length(candidate_substitutions{s}.labels);
            % approximate support of the pipeline's minimally required channels/types
            provides_labels(s) = isempty(fast_setdiff(lower(required_channels),lower(candidate_substitutions{s}.labels)));
            matches_types(s) = isequal(sort(lower(required_types)),sort(lower(candidate_substitutions{s}.types)));
            provides_types(s) = isempty(fast_setdiff(lower(required_types),lower(candidate_substitutions{s}.types)));
            provides_chancount(s) = length(candidate_substitutions{s}.labels) >= length(required_channels);
            % calculate overall score
            match_score(s) = matches_labels_and_order(s)*64 + matches_labels(s)*32 + provides_labels(s)*16 + matches_chancount(s)*8 + matches_types(s)*4 + provides_types(s)*2 + provides_chancount(s);
        end

        % decide which stream to use for this rawdata node
        [sorted_score,ranking] = sort(match_score,'descend'); %#ok<ASGLU>
        s = ranking(1); use_stream = candidate_substitutions{s};
        
        % inform about issues with the chosen match
        if ~matches_labels_and_order(s)
            if match_score(s) >= 16                
                message_header = 'NOTE: no perfectly matching stream found';
            else
                message_header = 'WARNING: only poorly matching streams found';
            end
            fprintf('%s; matches ordered labels: %i, matches unordered labels: %i, provides necessary labels: %i, matches channel count: %i, matches channel types: %i, provides necessary types: %i, provides necessary channel count: %i\n',...
                message_header,matches_labels_and_order(s),matches_labels(s),provides_labels(s),matches_chancount(s),matches_types(s),provides_types(s),provides_chancount(s));
            if ~matches_labels(s)
                fprintf('  pipeline input labels: %s\n  using stream''s labels: %s\n',hlp_tostring(use_stream.labels),hlp_tostring(pip.parts{1})); end
        end
        
        % warn about potential ambiguities
        if length(ranking) > 1
            same_rating = [false ranking(2:end)==ranking(1)];
            if any(same_rating)
                fprintf('WARNING: more than one stream matched the given pipeline equally well: %s\n', hlp_tostring(candidate_substitutions(same_rating),1000));
                second_candidate = candidate_substitutions(find(same_rating,1));
                fprintf('  other stream''s labels: %s\n',hlp_tostring(second_candidate.labels));
            end
        end

        if ischar(use_stream.substitute)
            % substitute a variable name into the rawdata node
            pip = struct('head',@rawdata,'parts',{{use_stream.substitute,1:length(use_stream.labels)}});
        else
            % substitute the rawdata node by the stream's expression
            pip = use_stream.substitute;
        end

        return; % done with substitution, no need to further recurse
        
    case 'flt_selchans'
        % a channel selection node: allows us to potentially determine the set of channel required from downstream nodes
        if isempty(required_channels)
            try
                args = hlp_microcache('resolve_streams',@arg_report,'vals',pip.head,pip.parts);
                if iscellstr(args.channels) && ~args.find_closest && ~args.remove_selection && ~args.relabel_to_query
                    required_channels = args.channels; end
                if iscellstr(args.channels) && (args.find_closest  || args.relabel_to_query)
                    required_channels = {}; end
            catch e
                disp_once('Error while parsing a flt_selchans node (syntax changed?): %s',e.message);
            end
        end
        
    case 'flt_seltypes'
        % a type selection node: allows us to potentially determine the set of types required from downstream nodes
        if isempty(required_types)
            try
                args = hlp_microcache('resolve_streams',@arg_report,'vals',pip.head,pip.parts);
                if iscellstr(args.chantypes)
                    required_types = args.chantypes; end
            catch e
                disp_once('Error while parsing a flt_seltypes node (syntax changed?): %s',e.message);
            end
        end
        
    otherwise
        % a generic node: if we can be sure that this node treats channels independently, we can
        % propagate the required set of channels and types further downstream; in all other cases
        % we assume that the node requires all channels that the lower levels provide
        try
            props = hlp_microcache('props',@arg_report,'properties',pip.head,{hlp_fileinfo(pip.head,[],'hash')});
        catch e
            disp_once('The properties of filter stage "%s" could not be queried due to error: %s',char(pip.head),hlp_handleerror(e));
            props = [];
        end
        if ~isfield(props,'independent_channels') || ~props.independent_channels
            % the node doesn't have independent channels: need all channnels & types from below
            required_channels = {};
            required_types = {};
        end
end

% recurse into expression-typed subnodes
for p = find(cellfun(@(p)all(isfield(p,{'head','parts'})),pip.parts))
    pip.parts{p} = resolve_rawdata(pip.parts{p},candidate_substitutions,required_channels,required_types); end
