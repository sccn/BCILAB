function res = utl_partition_bundle(bundle,inds,varargin)
% Internal. Cross-validation partitioner for stream bundles.
%
% In:
%   Bundle   : a stream bundle
%
%   IndexSet : partitioner index set -- see set_partition
%
%   EpochBunds : optional upper bound on epochs -- see set_partition
%
% Out:
%   Result : result of the operation (either a partitioned bundle or an index set cardinality)
%
% See also:
%   set_partition, utl_default_partitioner
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-28
dp;

% input validation
if ~isstruct(bundle) || ~isscalar(bundle)
    error('The given Bundle argument must be a 1x1 struct, but was: %s',hlp_tostring(bundle,10000)); end
if ~isfield(bundle,'streams')
    error('The given Bundle argument must be a 1x1 struct, but was: %s',hlp_tostring(bundle,10000)); end
if ~iscell(bundle.streams) || isempty(bundle.streams)
    error('The given Bundle argument has an invalid .streams field (must be a nonempty cell array), but was %s',hlp_tostring(bundle.streams,10000)); end

if isempty(inds)
    % compute index set size (from first stream)
    res = exp_eval_optimized(set_partition(bundle.streams{1},[],varargin{:}));    
elseif isnumeric(inds)
    % validate indices
    if min(size(inds)) ~= 1
        error('The given IndexSet argument must be a vector, but was: %s',hlp_tostring(inds,10000)); end
    if ~(all(round(inds) == inds) && all(inds > 0))
        error('The given index vector must contain positive integers, but was: %s',hlp_tostring(inds,10000)); end
    % partition the streams (symbolically)
    for s=1:length(bundle.streams)
        bundle.streams{s} = set_partition(bundle.streams{s},inds,varargin{:}); end    
    res = bundle;
else
    error('The given IndexSet argument is malformed: must be a vector of indices, but was: %s',hlp_tostring(inds,10000));
end
