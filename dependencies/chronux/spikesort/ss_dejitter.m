function spikes = ss_dejitter(spikes, maxshift)

% SS_DEJITTER  Aligns waveform peaks.
%     SPIKES = SS_DEJITTER(SPIKES) takes and returns a spike-sorting object
%     SPIKES and dejitters the spike waveforms using a 'center-of-mass'
%     method with a maximum shift of 3 (see below).
%
%     'Dejittering' refers to the registration of digitally sampled/thresholded 
%     waveforms.  The sampling/thresholding procedure introduces noise due to
%     (1) variation in relative timing between an analog waveform and the
%     sample clock and (2) variation in the time of threshold crossing due to
%     noise.  The resulting misalignment even in otherwise identical analog 
%     waveforms is known as jitter.  Dejittering involves estimating the location
%     of a fiducial (e.g., the central peak) for each waveform and aligning these
%     markers across waveforms.
%
%     SPIKES = SS_DEJITTER(SPIKES, MAXSHIFT) also takes MAXSHIFT argument
%     (default: 3), specifying the maximum shift allowed during dejittering.
% 
%     The number of samples per waveform will be altered by this function when 
%     samples are made invalid by realignment (i.e., the amount of shift requires
%     extrapolating past the boundaries of the original waveform) and the returned
%     SPIKES object will have its waveforms clipped to exclude invalid regions.
%     The MAXSHIFT argument limits this by specifying the maximum allowed shift
%     (in samples).  Waveforms requiring more than a MAXSHIFT samples realignment
%     are taken to be outliers and will not be shifted.
% 
%     Dejittering requires the definition of threshold crossing time (threshT) and
%     threshold level vector (threshV) with both low and high thresholds (use +/-
%     Inf for either the low or the high threshold if it was not used in extracting
%     the waveforms).  The 'center-of-mass' method considers the region following
%     threshT that remains below/above threshold and computes
%          (sum_t[t * (v_t - thresh)] / sum_t[(v_t - thresh)]).
%     to obtain the location of the peak.  The interpolation to align the waveforms
%     is then performed using a cubic spline.
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%     Sahani M (1999).  PhD thesis, Pasadena, CA: Caltech.

%   Last Modified By: sbm on Thu Oct  6 20:30:26 2005

starttime = clock;

%%%%%%%%%% DEFAULTS
if (nargin < 2),  maxshift = 3;    end;   % default max shift

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
elseif (~isfield(spikes, 'threshT'))
	error('SS:threshT_undefined', 'The SS object must define the sample index of threshold crossing.');
elseif (spikes.threshT > size(spikes.waveforms, 2))
    error('SS:threshT_invalid', 'The threshold time ThreshT does not fall within the supplied data waveform.');
elseif (~isfield(spikes, 'threshV') || (length(spikes.threshV) ~= 2))
	error('SS:threshV_invalid', 'The SS object must define both low and high thresholds.  Use +/- Inf if only one of the two was used.'); 
elseif ((spikes.threshV(2) - spikes.threshV(1)) < 0)
    error('SS:threshV_illegal', 'The SS object high threshold must be greater than its low threshold.');
end


%%%%%%%%%% CONSTANTS
% spline_chunk = 5000;   % # of splines to do at a time; keeps memory use down
numspikes = size(spikes.waveforms, 1);
numsamples = size(spikes.waveforms, 2);

%%%%%%%%%% FIND FIDUCIALS
% first, get a mask over the peaks
isThreshLo = (spikes.waveforms(:, spikes.threshT) < spikes.threshV(1));
[h, w, pl, mask] = thresholded_peaks(spikes);

% next, subtract appropriate thresholds, making all values in the peak
% have positive sign, since we're about to treat them as 'mass'
waves = spikes.waveforms;
waves(~isThreshLo, :) = waves(~isThreshLo, :) - spikes.threshV(2);
waves( isThreshLo, :) = spikes.threshV(1) - waves( isThreshLo, :);
waves = waves .* mask;

% now COM is straightforward
fiducials = (waves * (1:numsamples)') ./ sum(waves, 2);

clear mask waves;  % manual garbage collection

%%%%%%%%%% REALIGN FIDUCIALS
% We line up the fiducials around the sample that requires the least shifting.
target = round(mean(fiducials));

% Determine shifted indices for each waveform
shifts = fiducials - target;
shifts(abs(shifts) > maxshift) = 0;   % big shifts are outliers; don't alter these
resample_inds = repmat((1:numsamples), [numspikes, 1]) + repmat(shifts, [1, numsamples]);

% Which regions are invalid?
left_valid = max(find(any(resample_inds < 1))) + 1;            if (isempty(left_valid)), left_valid = 1;  end;
right_valid = min(find(any(resample_inds > numsamples))) - 1;  if (isempty(right_valid)), right_valid = numsamples; end;
if (left_valid >= right_valid)
    error('SS:waveform_invalidated', 'Realignment invalidates all samples.  Try reducing MAXSHIFT.');
end

% Although 'spline' can find independent interpolants for each row of a matrix
% (warning: the help for 'spline' is misleading on rows vs. cols), it does not
% allow each row to be resampled on its own grid.  To get around this (since the
% for-loop approach is ~4x slower), we obtain the piecewise polynomial coefficients
% along each waveform, pad to handle endpoints, and then string them all together
% into one long series.  This will let us compute the shifted waveforms in one 
% fell swoop . . . the only caveat is that requests to extrapolate beyond
% [1:numsamples] become meaningless; but we're not interested in these anyway.
pp = spline((1:numsamples), spikes.waveforms);

% The coefficients from 'spline' come out ordered by column rather than by
% row/waveform, so we reorder.
pp.coefs = reshape(pp.coefs, numspikes, numsamples-1, []);
pp.coefs = permute(pp.coefs, [2 1 3]);

% DEBUGGING -- see the splie
% example = 1;
% figure;  cla;  plot(spikes.waveforms(example,:), 'm.-', 'LineWidth', 4);  hold on;
% for t0 = 1:(numsamples-1)
%     t = linspace(0,1,10);   v = polyval(squeeze(pp.coefs(t0,example,:)), t);  plot(t+t0,v,'k');
% end
% hold off;

% 'ppeval' uses the right spline at endpoints; i.e., if we just strung together
% the splines as is, the last sample of each waveform would be interpolated
% using the first polynomial from the next waveform.  So we add a DC polynomial
% to the end of each spike whose value equals to the last sample of the spike.
% Now after stringing together, interpolation at ((r-1)*numsamples + [1:numsamples])
% gives the rth spike (within roundoff error).
padzeros = zeros(1, numspikes, 4);
pp.coefs = cat(1, pp.coefs, padzeros);
pp.coefs(numsamples,:,4) = spikes.waveforms(:,end)';

% Make the 'strung-together' version.
pp.coefs = reshape(pp.coefs, [], 4);

% Rework the remaining information in the piecewise polynomial to be consistent.
pp.pieces = numsamples * numspikes;   pp.dim = 1;   pp.breaks = (1:(pp.pieces+1));

% Change resample indices to correspond to the strung-together piecewise polynomial.
offset = ((1:numspikes) - 1)' * numsamples;
resample_inds = resample_inds + repmat(offset, [1, numsamples]);

% OK, finally.  Do the resample & clip the invalid regions
resampled = ppval(pp, resample_inds);
spikes.waveforms = resampled(:, left_valid:right_valid);

% Also, the threshold crossing time has changed due to data clipping
spikes.threshT = spikes.threshT - left_valid + 1;

% Worse, the threshold time/value may have changed while shifting.  This
% will make it hard to identify peaks later, so we need new values.
meanlevel = mean(spikes.waveforms(:));   % estimate of baseline
notallowed = [1:(spikes.threshT-maxshift-1),(spikes.threshT+maxshift+1):(size(spikes.waveforms, 2))];

deflection = mean(abs(spikes.waveforms - meanlevel), 1);  % mean deviation from baseline
deflection(notallowed) = 0;
[junk, spikes.threshT] = max(deflection);  % max deflect near old threshT marks new threshold marker
threshValues = spikes.waveforms(:, spikes.threshT);  % values at new marker
spikes.threshV(1) = max([threshValues(threshValues < meanlevel); -Inf]); % value closest to baseline but below it
spikes.threshV(2) = min([threshValues(threshValues > meanlevel); Inf]); % value closest to baseline but above it 

spikes.tictoc.dejitter = etime(clock, starttime);
