function [allow, scores, cdfs] = isiQuality(unit1times, unit2times, tmin, tmax, tref, Fs)

% ISIQUALITY   Computes statistical measures of refactory period quality.  
%     [allow, scores, cdfs] = isiQuality(unit1times, unit2times, tmin, tmax, tref, Fs)
%
% Returns a boolean specifying whether the refractory period statistics are
%   statistically worsened by combining the two input lists, along with two
%   outputs that give more details:
%
% INPUTS:
%   unit1times, unit2times  : (sec) original spike time lists
%   tmin                    : (sec) minimum possible interval
%   tmax                    : (sec) max interval for statistical comparison 
%   tref                    : (sec) estimate of the refractory period
%   Fs                      : (Hz)  data sampling frequency
%
% OUTPUTS:
%   allow                   : (0 or 1) is it ok to combine the two lists?
%   scores                  : (3 x 1) isi scores for [unit1 unit2 combined] lists;
%                             describe the normalized fraction of spikes violating
%                             the putative refractory period (1.0 means that the
%                             mean density of intervals below tref is the same as the
%                             mean density from tref to tmax.
%   cdfs                    : (3 x m) isi cdfs for unit1, unit2 and combined lists,
%                             covering the range from 0 to tmax with bin size
%                             1/Fs. (i.e., m is floor(tmax*Fs)).

% Useful numbers
f = (tref - tmin) ./ (tmax - tmin);        % fraction of range that is below tref
bins = linspace(0, tmax, floor(tmax*Fs));  % using tmax*Fs bins means no info is lost to rounding
refractory_bin = sum((bins <= tref));      % index of refractory period

% Convert from times to interspike intervals
isi1 = diff(sort(unit1times));
isi2 = diff(sort(unit2times));
isiT = diff(sort([unit1times; unit2times]));

% Calculate histograms for the intervals below tmax.
isi1hist = hist(isi1(isi1 < tmax), bins);
isi2hist = hist(isi2(isi2 < tmax), bins);
isiThist = hist(isiT(isiT < tmax), bins);

% Count total # intervals below tmax (set 0 counts to 1 as a courtesy to the next step)
isi1count = max(sum(isi1hist),1);
isi2count = max(sum(isi2hist),1);
isiTcount = max(sum(isiThist),1);

% Make cdfs from histograms
cdfs = zeros(3, length(bins));
cdfs(1,:) = cumsum(isi1hist) ./ isi1count;
cdfs(2,:) = cumsum(isi2hist) ./ isi2count;
cdfs(3,:) = cumsum(isiThist) ./ isiTcount;

% 
% if ((isi1count == 1) || (isi2count == 1)),  % if either original list had no violations
% 
% 	
% else
	% Compute the (scaled) difference between the initial cdfs and the combined cdfs.
	diffs = zeros(2, length(bins));
	diffs(1,:) = sqrt((isi1count*isiTcount)./(isi1count+isiTcount)) .* (cdfs(3,:) - cdfs(1,:));
	diffs(2,:) = sqrt((isi2count*isiTcount)./(isi2count+isiTcount)) .* (cdfs(3,:) - cdfs(2,:));
	
	% Max value of this difference in the region shorter than the refractory period

	dstats(1) = max(diffs(1, 1:refractory_bin));
	dstats(2) = max(diffs(2, 1:refractory_bin));
	
	% based on Kolmogorov-Smirnoff statistics (see Fee, Mitra, Kleinfeld, 1996)
	lambda = 0:0.01:2;
	pdfDstat = (0.5 * (1 + erf(lambda./(sqrt(2*f*(1-f)))))) - ...
		(0.5 .* exp(-2.*lambda.^2) .* (1 - erf((1-2*f).*lambda./(sqrt(2*f*(1-f))))));
	cutoff = lambda(max(find(pdfDstat < 0.95)));

	% Does either d-statistic exceed the statistical cutoff?
	if (any(dstats > cutoff)),  allow = 0;
    else                        allow = 1;
	end
% end

% The cdf at the tref bin is the fraction of spikes below tref; to get a score,
% we normalize by dividing by f from above.
scores =  (cdfs(:, refractory_bin)) ./ f;
