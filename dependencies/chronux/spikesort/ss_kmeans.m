function spikes = ss_kmeans(spikes, options)
% SS_KMEANS  K-means clustering.
%     SPIKES = SS_KMEANS(SPIKES) takes and returns a spike-sorting object SPIKES.
%     It k-means clusters the (M x N) matrix SPIKES.WAVEFORMS and stores the
%     resulting group assignments in SPIKES.OVERCLUSTER.ASSIGNS, the cluster
%     centroids in SPIKES.OVERCLUSTER.CENTROIDS, and the mean squared error in
%     SPIKES.OVERCLUSTER.MSE.  The W, B, and T matrices (within-cluster, between-
%     cluster, and total) covariance matrices are in SPIKES.OVERCLUSTER.W, etc.
%
%     K-means algorithm is an EM-like algorithm that finds K cluster centers in
%     N-dimensional space and assigns each of the M data vectors to one of these
%     K points.  The process is iterative: new cluster centers are calculated as
%     the mean of their previously assigned vectors and vectors are then each 
%     reassigned to their closest cluster center.  This finds a local minimum for
%     the mean squared distance from each point to its cluster center; this is
%     also a (local) MLE for the model of K isotropic gaussian clusters.
% 
%     This method is highly outlier sensitive; MLE is not robust to the addition
%     of a few waveforms that are not like the others and will shift the 'true'
%     cluster centers (or add new ones) to account for these points.  See
%     SS_OUTLIERS for one solution.
% 
%     The algorithm used here speeds convergence by first solving for 2 means,
%     then using these means (slightly jittered) as starting points for a 4-means
%     solution.  This continues for log2(K) steps until K-means have been found.
%
%     SPIKES = SS_KMEANS(SPIKES, OPTIONS) allows specification of clustering
%     parameters.  OPTIONS is a structure with some/all of the following fields
%     defined.  (Any OPTIONS fields left undefined (or all fields if no OPTIONS
%     structure is passed in) uses its default value.)
%
%         OPTIONS.DIVISIONS (default: round(log2(M/500)), clipped to be between
%               4 and 8) sets the desired number of clusters to 2^DIVISIONS.  The
%               actual number of clusters may be slightly more/less than this #.
%         OPTIONS.REPS (default: 1) specifies the number of runs of the full
%               k-means solution.  The function will return the assignments that
%               resulted in the minimum MSE.
%         OPTIONS.REASSIGN_CONVERGE (default: 0) defines a convergence condition
%               by specifying the max number of vectors allowed to be reassigned
%               in an EM step.  If <= this number of vectors is reassigned, the
%               this condition is met.
%         OPTIONS.MSE_CONVERGE (default: 0) defines a second convergence condition.
%               If the fractional change in mean squared error from one iteration
%               to the next is smaller than this value, this condition is met.
%
%     NOTE: Iteration stops when either of the convergence conditions is met.
%  
% References:
%     Duda RO et al (2001).  _Pattern Classification_, Wiley-Interscience

%   Last Modified By: sbm on Sun Aug 13 02:28:36 2006

% Undocumented options: 
%       OPTIONS.PROGRESS (default: 1) determines whether the progress
%                                     bar is displayed during the clustering.
%       OPTIONS.SPLITMORE (default: 1) determines whether to do extra cluster
%                                       splitting to try for smaller clusters

starttime = clock;

% save the random number seeds before we do anything, in case we want to
% replicate this run exactly ...
randnstate = randn('state');  randstate = rand('state');

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
end

%%%%%%%%%% CONSTANTS
target_clustersize = 500;
waves = spikes.waveforms;         % ref w/o struct is better for R13 acceleration
[M,N] = size(waves);
jitter = meandist_estim(waves) / 100 / N;        % heuristic

%%%%%%%%%% DEFAULTS
opts.divisions = round(log2(M / target_clustersize));  % power of 2 that gives closest to target_clustersize, but in [4..7]
opts.divisions = max(min(opts.divisions, 8), 4);       %     restrict to 16-256 clusters; heuristic.
opts.memorylimit = M;                   % all spikes at once
opts.memorylimit = 300;                 % 
opts.reps = 1;                          % just one repetition
opts.reassign_converge = 0;             % the last few stubborn ones are probably borderline cases anyway ...
opts.reassign_rough = round(0.005*M);   % (no need at all to get full convergence for intermediate runs) 
opts.mse_converge = 0;                  % ... and by default, we don't use the mse convergence criterion
opts.progress = 1;
opts.splitmore = 1;
if (nargin > 1)
	supplied = lower(fieldnames(options));   % which options did the user specify?
	for op = 1:length(supplied)              % copy those over the defaults
		opts.(supplied{op}) = options.(supplied{op});  % this is the preferred syntax as of R13 --
	end
end

%%%%%%%%%% CLUSTERING
normsq = sum(waves.^2, 2);
assigns = ones(M, opts.reps);
mse = Inf * ones(1, opts.reps);
clear fast_kmeans;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT NOTE: For BLAS efficiency reasons,  %
%   the waveforms & centroids matrices are      %
%   transposed from their ordinary orientation  %
%   in the following section of code.           %
%   E.g., waveforms is (D samples x N waveforms)%
waves = waves';                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for rep = 1:opts.reps                                 % TOTAL # REPETITIONS
	
	centroid = mean(waves, 2);  % always start here
    itercounter = zeros(1,opts.divisions);
	for iter = 1:opts.divisions                       % # K-MEANS SPLITS
		itercounter(iter) = 0;
		oldassigns = zeros(M, 1);  oldmse = Inf;
		assign_converge = 0;       mse_converge = 0;
        
        if (iter == opts.divisions)
            progress = opts.progress;  reassign_criterion = opts.reassign_converge; 
        else
            progress = 0;              reassign_criterion = opts.reassign_rough;
        end

        if ((iter==opts.divisions) && opts.splitmore)  % do some extra splitting on clusters that look too big
            toobig = find(clustersizes > 2*target_clustersize);
            splitbig = centroid(:,toobig) + randn(length(toobig),size(centroid,1))';
            centroid = [centroid splitbig];
        end

		centroid = [centroid centroid] + jitter * randn(2*size(centroid, 2),size(centroid,1))'; % split & jitter
		clustersizes = ones(1,size(centroid,2));  % temporary placeholder to indicate no empty clusters at first
        
        if (opts.progress)
            progressBar(0, 1, sprintf('Calculating %d means.', size(centroid,2))); % crude ...
        end

		while (~(assign_converge || mse_converge))     % convergence?
			%%%%% Clean out empty clusters
			centroid(:,(clustersizes == 0)) = [];
			            
			%%%%% EM STEP
			[assignlist, bestdists, clustersizes, centroid] = ...
				        fast_kmeans_step(waves, centroid, normsq, opts.memorylimit);
						
			%%%%% Compute convergence info
			mse(rep) = mean(bestdists);
            mse_converge = ((1 - (mse(rep)/oldmse)) <= opts.mse_converge);   % fractional change
			oldmse = mse(rep);
            
			changed_assigns = sum(assignlist ~= oldassigns);
			assign_converge = (changed_assigns <= reassign_criterion);   % num waveforms reassigned
            if (progress)
                progressBar(((M - changed_assigns)/M).^10, 5); % rough ...
            end
			oldassigns = assignlist;
			itercounter(iter) = itercounter(iter) + 1;
		end
	end
	
	% Finally, reassign spikes in singleton clusters to next best ...
	junkclust = find((clustersizes <= 1));   junkspikes = ismember(assignlist, junkclust);
	centroid(:,junkclust) = Inf;  % take these out of the running & recompute distances
    if (any(junkspikes))
        dist_to_rest = pairdist(waves(:,junkspikes)', centroid', 'nosqrt');
        [bestdists(junkspikes), assignlist(junkspikes)] = min(dist_to_rest, [], 2);
    end
	mse(rep) = mean(bestdists);
	
	assigns(:,rep) = assignlist;
end

spikes.overcluster.iteration_count = itercounter;
spikes.overcluster.randn_state = randnstate;
spikes.overcluster.rand_state = randstate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waves & centroids now return to their usual   %
% orientation.  E.g., waveforms is again        %
%    N waveforms x D samples                    %
waves = waves';                                 %
centroid = centroid';                           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Finish up by selecting the lowest mse over repetitions.
[bestmse, choice] = min(mse);
spikes.overcluster.assigns = sortAssignments(assigns(:,choice));
spikes.overcluster.mse = bestmse;
spikes.overcluster.sqerr = bestdists;

% We also save the winning cluster centers as a convenience
numclusts = max(spikes.overcluster.assigns);
spikes.overcluster.centroids = zeros(numclusts, N);
for clust = 1:numclusts
	members = find(spikes.overcluster.assigns == clust);
    spikes.overcluster.centroids(clust,:) = mean(waves(members,:), 1);
end

% And W, B, T matrices -- easy since T & B are fast to compute and T = W + B)
spikes.overcluster.T = cov(waves);             % normalize everything by M-1, not M
spikes.overcluster.B = cov(spikes.overcluster.centroids(spikes.overcluster.assigns, :));
spikes.overcluster.W = spikes.overcluster.T - spikes.overcluster.B;

% Finally, assign colors to the spike clusters here for consistency ...
cmap = jetm(numclusts);
spikes.overcluster.colors = cmap(randperm(numclusts),:);

if (opts.progress), progressBar(1.0, 1, ''); end
spikes.tictoc.kmeans = etime(clock, starttime);
