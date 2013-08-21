function spikes = merge_clusters(spikes, to, from)
% MERGE_CLUSTERS  Merge two clusters after automatic hierarchical aggregation.
%    SPIKES = MERGE_CLUSTERS(SPIKES, TO, FROM) takes and returns a spike-
%    sorting object SPIKES.  SPIKES must have gone through a hierarchical
%    clustering aggregation (e.g., SS_AGGREGATE) previous to this function call.
%   
%    All spikes belonging to the cluster whose label number is given by FROM
%    are merged into the cluster with label TO.  The hierarchical clustering
%    assignments and aggregation tree are modified to reflect this change.
%    If SPIKES contains a SPIKETIMES vector, the ratio of the interspike  
%    interval count below 2 msec to the count below 10 msec is computed and
%    entered into the SPIKES.HIERARCHY.TREE.

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'hierarchy'))
	error('SS:hierarchical_information_unavailable', 'Hierarchical clustering must be performed before attempting to merge clusters.');
elseif (~all(ismember([to,from], unique(spikes.hierarchy.assigns))))
    error('SS:cluster_numbers_not_found', 'One or both of the cluster labels supplied does not exist.');
elseif ((length(from) > 1) || (length(to) > 1))
    error('SS:one_at_a_time', 'The ''from'' and ''to'' labels must be scalar values.');
end
if (to == 0)
    warning('SS:merge_outliers', 'Adding spikes to the outlier cluster.');
end

%%%%%%%%%% MERGING ASSIGNMENTS
members_from = find(spikes.hierarchy.assigns == from);    % Get list of spikes to move ...
orig_members_to = find(spikes.hierarchy.assigns == to);   %   (we need this below)
spikes.hierarchy.assigns(members_from) = to;              % ... and relabel them.

%%%%%%%%%% COMPUTE CONNECTION STRENGTHS FROM INTERFACE ENERGY
%% Temporary hack; this should be recomputed from the interface energy
cs = 1;

%%%%%%%%%% MODIFY AGGREGATION TREE
tmin = size(spikes.waveforms,2)./spikes.Fs;  tref = max(0.002, tmin*1.5);
if (isfield(spikes, 'spiketimes'))   % might as well compute the isi score if we can ...
    t1 = spikes.spiketimes(orig_members_to);
    t2 = spikes.spiketimes(members_from);
    [allow, score] = isiQuality(t1, t2, tmin, 0.010, tref, spikes.Fs);
    score = score(3);
else                                 % ... but no stress if we can't
    score = 0;
end
spikes.hierarchy.tree = [spikes.hierarchy.tree; to from cs score];

