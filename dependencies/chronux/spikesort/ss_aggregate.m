function spikes = ss_aggregate(spikes, reintegrate_outliers)
% SS_AGGREGATE  ISI-restricted heirarchical cluster aggregation.
%     SPIKES = SS_AGGREGATE(SPIKES) takes and returns a spike-sorting object
%     SPIKES after aggregating clusters (requires a previous overclustering
%     and an interface energy matrix calculation).  The aggregation tree is
%     stored in SPIKES.HIERARCHY.TREE and the new assignments are stored in
%     SPIKES.HIERARCHY.ASSIGNS.
%
%     The algorithm computes a similarity/connection matrix using the interface
%     energy.  It then chooses the cluster pair with the highest connection
%     strength, aggregates them, recalculates connection strengths, and then
%     repeats the process.  Cluster aggregation is contingent on passing an
%     interspike interval (ISI) test if SPIKES.SPIKETIMES is defined; if this
%     test is not passed, the pair is not aggregated and aggregation continues.
%     Aggregation stops when either (1) all remaining pairs fail the ISI test
%     or (2) the connection strength drops below a (heuristic) cutoff of 0.01.
%
%     The ISI test needs an idea of the expected refractory period to test
%     for violations; as a default, it uses 1.5 times the width of the
%     waveforms in the SPIKES structure.  To override this, define the
%     value SPIKES.OPTIONS.REFRACTORY_PERIOD with units of seconds (e.g.,
%     set to 0.0017 to indicate an expected refractory period of 1.7 msec).
%
%     The SPIKES.HIERARCHY.TREE output is a matrix describing the aggregation.
%     Each aggregation step entry produces a row, listed in the order they
%     were performed.  The first two columns are the indices of the clusters
%     that were aggregated; the index assigned to the aggregate for future
%     entries is the lower of the two original indices.  The third column is
%     the connection strength between the clusters before aggregation and the
%     fourth column is the isi statistic for the aggregate (0 if isi statistics
%     are not being used).
% 
%     After aggregation, outliers that were previously removed are typically 
%     reinserted into the spikes list so that the list aligns with the original
%     (pre-sorted) list.  The outlier waveforms and spike times are thus by
%     default added back to the SPIKES.WAVEFORMS and SPIKES.SPIKETIMES fields
%     respectively, after all other aggregation is complete.  These waveforms
%     are assigned the label 0 in SPIKES.HIERARCHY.ASSIGNS and the other, 
%     non-outlier, spikes are renumbered accordingly.  To prevent this from
%     occuring, pass in 0 as the second argument to this function, i.e.,
%     SPIKES = SS_AGGREGATE(SPIKES, 0);  
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88

debug = 1;

starttime = clock;

cutoff = 0.01;   % arbitrarily stop aggregation when overlap density is < 1% of main cluster density

s = warning('off', 'MATLAB:divideByZero');

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'hierarchy') || ~isfield(spikes.hierarchy, 'interface_energy'))
    error('SS:energy_not_computed', 'An energy matrix must be computed before aggregation.');
elseif (~isfield(spikes, 'overcluster') || ~isfield(spikes.overcluster, 'assigns'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before aggregation.');
elseif (isfield(spikes, 'spiketimes') && ~isfield(spikes, 'Fs'))
    error('SS:bad_stop_condition', 'A sampling frequency Fs must be supplied for an ISI stop condition.');
end
if (isfield(spikes.hierarchy, 'assigns') && any(spikes.hierarchy.assigns == 0))
    error('SS:aggregate_after_outliers', 'Aggregation can not be performed after outliers are reintegrated into the data.');
end
if (nargin < 2),  reintegrate_outliers = 1;  end

tmin = size(spikes.waveforms,2)./spikes.Fs;  % minimum possible time btw spikes
if (isfield(spikes, 'options') && isfield(spikes.options, 'refractory_period'))
    tref = spikes.options.refractory_period;
else
    tref = max(0.002, tmin*1.5);                 % crude guess as to refractory period
end

%%%%% INITIALIZE A FEW THINGS
assignments = spikes.overcluster.assigns;
interface_energy = spikes.hierarchy.interface_energy;
numclusts = max(assignments);
numpts = full(sparse(assignments, 1, 1, numclusts, 1));
tree = [];
untested = 3*ones(numclusts);    % they're all initially untested

% Energy merging progress grid ...
handle_fig = figure;  handle_img = imagesc(untested);  axis square;
colormap([0 0 0; 0.9 0.4 0.4; 1 1 0; 0.4 0.8 0.2]);  % [bk, rd, yl, gn] => (0 combined, 1 not allowed, 2 testing, 3 untested)

%%%%% AGGREGATE HIGHEST CONNECTION STRENGTHS UNTIL ALL TRIED
while (any(any(triu(untested,1)==3)))   % only the upper triangular half is meaningful
    % compute connection strengths from interface energies
    %   first, normalize energy:
    normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
    normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
    norm_energy = interface_energy ./ normalize;
    %   then, compute connection strength
    self = repmat(diag(norm_energy), [1,numclusts]);
    connect_strength = 2 .* norm_energy ./ (self + self');
    connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves

    % Find best remaining pair
    remaining = ((untested == 3) .* connect_strength);
    best = max(remaining(:));           % highest untested connection strength
    
    if (best < cutoff)   % No point continuing if connection strengths have gotten really lousy
        break;
    end
    
    [clust1 clust2] = find(connect_strength == best);  % who're the lucky winners?
    first = min(clust1(1),clust2(1));   % if we get 2 best pairs, just take the 1st
    second = max(clust1(1),clust2(1)); 
    untested(first,second) = 2;    untested(first,second) = 2;     % mark that we're trying this one
    set(handle_img, 'CData', triu(untested)); title(['Trying ' num2str(first) ' and ' num2str(second)]); drawnow;
	
    % Is this aggregation allowed?
    if (isfield(spikes, 'spiketimes'))  % if we were given spike times, use them ...
        t1 = spikes.spiketimes(assignments == first);
        t2 = spikes.spiketimes(assignments == second);    
        if (debug)
            score = sum((diff(sort([t1; t2])) < tref)) ./ (length(t1) + length(t2));
            allow = (score < 0.05);
            scores = [0 0 score];
        else
            [allow, scores] = isiQuality(t1, t2, tmin, 0.010, tref, spikes.Fs);
        end        
    else  % ... otherwise, there are no restrictions on aggregation
        allow = 1;
        scores = [0 0 0];
    end

    if (allow)      % Bookkeeping ...
        % Aggregation subsumes the higher index cluster into the lower.  Start by adding
        % (denormalized) interaction energies for the second (higher index) cluster
        % to those of the first and zeroing the old entries of the second.  Because of the
        % triangular structure of the connection matrix, repeat for both rows and columns ...
        interface_energy(first,:) = interface_energy(first,:) + interface_energy(second,:);
        interface_energy(second,:) = 0;
        interface_energy(:,first) = interface_energy(:,first) + interface_energy(:,second);
        interface_energy(:,second) = 0;
        interface_energy(second,second) = 1;  % keep self-energy at 1 (we may divide by it later)
        % since we added rows & columns, some energy values will have spilled over into the
        % lower half of the energy matrix (which must be upper triangular).  The next 2 steps
        % recover those values.
        overflow = tril(interface_energy, -1);   % spillover below diagonal
        interface_energy = interface_energy + overflow' - overflow;  % reflect above diagonal

        % update counts vector
        numpts(first) = numpts(first) + numpts(second);
        numpts(second) = 2;   % leaving this as 2 prevents div by zero during normalization above
        
        % make a tree entry for the aggregation we just performed
        tree = cat(1, tree, [first, second, best, scores(3)]);

        % Now actually change the numbers
        assignments(assignments == second) = first;
        
        % Finally, indicate that potential aggregations between the new cluster and 
        % other (nonempty) clusters are untested while pairs involving clusters that
        % have already been emptied should not be tested.
        untested(first,:) = 3;               untested(:,first) = 3;
        untested(tree(:,2),:) = 0;           untested(:,tree(:,2)) = 0;
    else
        untested(first,second) = 1;          untested(second,first) = 1;
    end
end
close(handle_fig);

%spikes.hierarchy.interface_energy_aggregated = interface_energy;
spikes.hierarchy.tree = tree;
spikes.hierarchy.assigns = assignments;

if (reintegrate_outliers && isfield(spikes, 'outliers') && ~isempty(spikes.outliers.badinds))
    % First, we make room by putting all of the non-outliers back into their original places
    spikes.waveforms(spikes.outliers.goodinds,:) = spikes.waveforms;
    spikes.spiketimes(spikes.outliers.goodinds,:) = spikes.spiketimes;
    spikes.hierarchy.assigns(spikes.outliers.goodinds) = spikes.hierarchy.assigns;
    
    % Then we fill in the outliers ...
    spikes.waveforms(spikes.outliers.badinds,:) = spikes.outliers.waveforms;
    spikes.spiketimes(spikes.outliers.badinds,:) = spikes.outliers.spiketimes;
    spikes.hierarchy.assigns(spikes.outliers.badinds) = 0;  % ... and add the '0' label.
    
    % We'll also want to add the assignments to the 'overcluster' list (this is
    % important for post-clustering splitting).
    spikes.overcluster.assigns(spikes.outliers.goodinds) = spikes.overcluster.assigns;
    spikes.overcluster.assigns(spikes.outliers.badinds) = 0;
    
    spikes = rmfield(spikes, 'outliers');  % don't need this any more -- its redundant.
end

spikes.tictoc.aggregate = etime(clock, starttime);


warning(s);