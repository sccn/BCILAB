function spikes = ss_energy(spikes)
% SS_ENERGY  Interface energy based cluster similarity computation.
%     SPIKES = SS_ENERGY(SPIKES) adds an interface energy matrix to a
%     spike-sorting object in SPIKES.HIERARCHY.INTERFACE_ENERGY.
% 
%     The energy similarity matrix is calculated by applying an exponential
%     decay to all pairwise euclidean distances between waveforms from two
%     clusters (or within a single cluster for intra-cluster energy) and 
%     summing these distances.
%
%     The calculation ignores the energy due to the zero distance between
%     points and themselves; this removes a dependence of the density on
%     the absolute size of the cluster.  As a result, singleton clusters
%     do not have a well-defined energy and will cause an error.
%  
%     When each entry is normalized by the number of distinct contributing
%     pairs (Na*Nb for off diagonal entries and Na*(Na-1)/2 on the diagonal),
%     it approximates the fraction of pairs in a given cluster whose distance
%     is not much greater than the length constant of the exponential and thus
%     provides an estimate of local density.  This function does not, however,
%     normalize SPIKES.HIERARCHY.INTERFACE_ENERGY, since the normalized form is
%     inconvenient during cluster aggregation.  The normalization can readily
%     be done, however, with
%          normalize = ((numpts * numpts') - diag(numpts));
%          normalize = normalize - diag(0.5*diag(normalize));
%          normalized_energy = interface_energy ./ normalize;
%     where 'numpts' is a vector of cluster sizes.
%
%     The unnormalized energy matrix can be updated during aggregation without
%     the need to recompute it from scratch.  The intra-cluster energy E(AB,AB)
%     of a cluster AB formed by aggregating clusters A and B is given by
%              E(AB,AB) = E(A,A) + E(B,B) + E(A,B)
%     and the inter-cluster energy between any cluster C and an aggregate AB is
%                 E(AB,C) = E(A,C) + E(B,C)
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88

%   Last Modified By: sbm on Fri Oct  7 21:35:16 2005

starttime = clock;

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
elseif (~isfield(spikes, 'overcluster'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before computing energy');
end
numclusts = length(unique(spikes.overcluster.assigns));
waves = spikes.waveforms;

%%%%% PREPARE SOME INFORMATION
normsqr = sum(waves.^2,2);
pts = cell(numclusts,1);    % collect spike indices for each cluster
for clust = 1:numclusts   
    pts{clust} = find(spikes.overcluster.assigns == clust);
end
numpts = cellfun('length', pts);
if (any(numpts < 2))
    error('SS:energy_ill_defined', 'Clusters with fewer than 2 points do not have a defined energy.');
end

%%%%% HEURISTIC DISTANCE SCALE that seems to work.  The calculation is not too sensitive to this parameter.
scale = sqrt(sum(diag(spikes.overcluster.W))) / 10;

%%%%% PREPARE TO LOOP
total = (numclusts^2 + numclusts) / 2;
k = 1;
progressBar(0, max(floor(total/100),1), 'Computing Interaction Energies . . .')
interface_energy = zeros(numclusts);

%%%%% PAIRWISE DISTANCES LOOP
assigns = spikes.overcluster.assigns;
for clust1 = 1:numclusts
	% Deal with self-case first, because PAIRDIST works better with a
	% different syntax in this case.
	dists = pairdist(waves(pts{clust1},:), 'reuse');
	interface_energy(clust1,clust1) = fast_interface_energy(dists,scale);
	k = k + 1;
    for clust2 = (clust1+1):numclusts   % now for the rest ...
        dists = pairdist(waves(pts{clust1},:), waves(pts{clust2},:), 'reuse');
        interface_energy(clust1,clust2) = fast_interface_energy(dists,scale);
        k = k + 1;
        progressBar(k/total);
    end
end

%%%%% CORRECTION TERMS
% The energy matrix so far includes a contribution in the intra-cluster
% energies that is not found in the inter-cluster energies; namely, the
% computation of   sum_(all x) sum_(all y) e^(-dist/scale)   for
% intra-cluster energy includes cases where x == y (so dist == 0).
interface_energy = interface_energy - diag(numpts);     % So subtract this out.

% Also, we've double counted pairs in the intra-energy case, since dist(a,b)
% and dist(b,a) are not treated as distinct;
interface_energy = interface_energy - diag(0.5*diag(interface_energy));

%%%%% FINISH UP
spikes.hierarchy.interface_energy = interface_energy;
spikes.tictoc.energy = etime(clock, starttime);
