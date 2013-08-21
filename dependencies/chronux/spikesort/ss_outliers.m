function spikes = ss_outliers(spikes, reps)
% SS_OUTLIERS  K-means based outlier detection.
%     SPIKES = SS_OUTLIERS(SPIKES) takes and returns a spike-sorting object
%     SPIKES.  It identifies likely outliers in the SPIKES.WAVEFORMS and moves
%     them from the (initially M x N) SPIKES.WAVEFORMS matrix to a (P x N)
%     SPIKES.OUTLIERS.WAVEFORMS (usually P << M), removing them the original
%     SPIKES.WAVEFORMS matrix.  If the SPIKES.OUTLIERS.WAVEFORMS matrix already
%     exists, new outliers are added to the existing matrix.
%
%     Outlier spike timestamps are recorded in SPIKES.OUTLIERS.SPIKETIMES
%     (these event times are removed from SPIKES.SPIKETIMES).  The cell matrix
%     SPIKES.OUTLIERS.WHY stores a string describing the reasons that the
%     corresponding waveforms was marked as an outlier.
% 
%     The SS_OUTLIERS function identifies outliers using an ad hoc heuristic
%     based on a k-means clustering; waveforms that end up very far from their 
%     assigned centroid are likely outliers.   The scale for this determination
%     comes from the average of the covariance matrices of each cluster (i.e.,
%     1/(N-1) times the (N x N) within-group sum of squares matrix).  We take
%     this as an approximation to the noise covariance in outlier-free clusters 
%     and note that if the noise were locally Gaussian, then the waveform Mahalanobis
%     distances to assigned cluster means should be roughly Chi^2 distributed
%     (actually, F-distributed but we ignore the refinement for now).
%     NOTE: Outliers damage the k-means solution; the clustering is not robust to
%            gross violations of its (local) Gaussian assumption.  Repeat clustering  
%            on the cleaned data will thus yield a new solution ... which might
%            uncover further outliers.  This function attempts 3 cluster/clean
%            iterations (rule of thumb; tradeoff btw cleaning and run-time),
%            although it stops after any iteration that does not find an outlier.
%
%     The default Chi^2 CDF cutoff is (1 - 1/(100*M)).  To use a different cutoff,
%     specify the value in 'SPIKES.OPTIONS.OUTLIER_CUTOFF'; the closer the value
%     (e.g., 1-1e-8) is to 1, the fewer spikes will be considered outliers.
%
%     SPIKES = SS_OUTLIERS(SPIKES, REPS) performs REPS cluster/clean iterations
%     (default: 3).  
%
%     The reason given for these outliers in SPIKES.OUTLIERS.WHY is "K-Means Outliers".
%
%     NOTE: This function performs a k-means clustering and thus overwrites
%           existing k-means assignments in the SPIKES object.

%   Last Modified By: sbm on Thu Oct  6 20:15:56 2005

starttime = clock;

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
end

[M,N] = size(spikes.waveforms);
if (isfield(spikes, 'options') && isfield(spikes.options, 'outlier_cutoff'))
    cutoff = spikes.options.outlier_cutoff;
else
    cutoff = (1 - (1/(100*M)));    
end

if (nargin < 3),  reps = 3;  end
times_defined = isfield(spikes, 'spiketimes');

% Initialize the outliers sub-structure
if (~isfield(spikes, 'outliers'))
    spikes.outliers.waveforms = [];
    spikes.outliers.why = {};
    if (times_defined)
        spikes.outliers.spiketimes = [];
    end
    
    spikes.outliers.goodinds = (1:M)';  % We need these to re-insert the outlier 'cluster'
    spikes.outliers.badinds = [];      % into the waveforms matrix after sorting.
end

for cleaning = 1:reps   % Cluster/clean then rinse/repeat.  3 reps does a good job w/o taking too much time

    progressBar((cleaning-1)./reps, 1, ['Removing Outliers: Pass ' num2str(cleaning)]);
    
    %%%%% PERFORM A K-MEANS CLUSTERING OF THE DATA
    opts.mse_converge = 0.001;        % rough clustering is fine
    opts.progress = 0;
    spikes = ss_kmeans(spikes, opts);
    
    %%%%% MAHALANOBIS DISTANCES:   (x - x_mean)' * S^-1 * (x - x_mean)
    vectors_to_centers = spikes.waveforms - spikes.overcluster.centroids(spikes.overcluster.assigns,:);
    mahaldists = sum((vectors_to_centers .* (spikes.overcluster.W \ vectors_to_centers')'), 2)';
    
    %%%%% SPLIT OFF OUTLIERS
    bad = find(mahaldists > chi2inv(cutoff, N));   % find putative outliers
    
    if (isempty(bad))
        break;              % didn't find anything ... no sense continuing
    else
        %%%%% ADD OUTLIERS TO SS OBJECT AND REMOVE FROM MAIN LIST
        spikes.outliers.waveforms = cat(1, spikes.outliers.waveforms, spikes.waveforms(bad,:));
        spikes.outliers.why = cat(1, spikes.outliers.why, repmat({'K-means Outliers'}, [length(bad), 1]));
        spikes.outliers.badinds = cat(1, spikes.outliers.badinds, spikes.outliers.goodinds(bad(:)));
        if (times_defined)
            spikes.outliers.spiketimes = cat(1, spikes.outliers.spiketimes, spikes.spiketimes(bad,:));
            spikes.spiketimes(bad,:) = [];
        end
        spikes.outliers.goodinds(bad) = [];
        spikes.waveforms(bad,:) = [];
        spikes = rmfield(spikes, 'overcluster');  % this isn't supposed to be public-visible
    end
end

progressBar(1.0, 1, '');

%%%%% TIMING INFORMATION
spikes.tictoc = rmfield(spikes.tictoc, 'kmeans');
spikes.tictoc.outliers = etime(clock, starttime);
