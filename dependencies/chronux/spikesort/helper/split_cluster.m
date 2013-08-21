function spikes = split_cluster(spikes, clustnum)
% SPLIT_CLUSTERS  Split a cluster after automatic hierarchical aggregation.
%    SPIKES = SPLIT_CLUSTERS(SPIKES, CLUSTER_NUMBER) takes and returns a spike-
%    sorting object SPIKES.  SPIKES must have gone through a hierarchical
%    clustering aggregation (e.g., SS_AGGREGATE) previous to this function call.
%   
%    The spikes belonging to the cluster whose label number is given by
%    CLUSTER_NUMBER are split into the two clusters by undoing the last step
%    in the aggregation tree involving the cluster with label CLUSTER_NUMBER.
%    The hierarchical cluster assignments and aggregation tree is modified to
%    reflect this change.

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'hierarchy'))
	error('SS:hierarchical_information_unavailable', 'Hierarchical clustering must be performed before attempting to merge clusters.');
elseif (~ismember(clustnum, unique(spikes.hierarchy.assigns)))
    error('SS:cluster_number_not_found', 'The cluster label supplied does not exist.');
end

%%%%%%%%%% FINDING INDICES OF SPIKES TO SPLIT
tree = spikes.hierarchy.tree;                       % (convenient shorthand)
treeentry = max(find(tree(:,1) == clustnum));       % where was the target cluster last merged into ...
breakaway = tree(treeentry, 2);                     % ... and who was the lucky mergee?

% The breakaway cluster may have itself been created by merging several original clusters, so
% we need to walk up the aggregation list to find all original (i.e., overcluster) labels
% that contribute to newly formed breakaway.
for entry = (treeentry-1):-1:1
    if (ismember(tree(entry, 1), breakaway))
        breakaway = [breakaway tree(entry,2)];
    end
end

% Now get the list of spikes indices for the new cluster
members_breakaway = find(ismember(spikes.overcluster.assigns, breakaway));

%%%%%%%%%% DO THE SPLIT
spikes.hierarchy.assigns(members_breakaway) = breakaway(1);  % Using this label keeps the tree accurate
spikes.hierarchy.tree(treeentry,:) = [];
