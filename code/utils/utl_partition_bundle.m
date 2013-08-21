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

if isempty(inds)
    % compute index set size (from first stream)
    res = exp_eval(set_partition(bundle.streams{1},[],varargin{:}));
else
    % partition the streams (symbolically)
    for s=1:length(bundle.streams)
        bundle.streams{s} = set_partition(bundle.streams{s},inds,varargin{:}); end
    res = bundle;
end
