function [weights, anglebins] = anglehist(Z, bins)
%ANGLEHIST         Magnitude-weighted angle histogram of polar data.
%   WEIGHTS = ANGLEHIST(Z) histograms the angles of complex data Z into 10
%   equally spaced bins around 2*pi, with each element in Z contributing
%   an amount equal to its magnitude to its corresponding bin.
%
%   WEIGHTS = ANGLEHIST(Z, BINS), where BINS is a scalar, uses BINS bins.
%
%   [WEIGHTS, ANGLEBINS] = ANGLEHIST(...) returns the bin centers
%   corresponding to the weights.
%
%   ANGLEHIST(...) without output arguments plots a polar bar histogram
%   of the results.

% Convert data into a column vector
Z = Z(:);

% Argument parsing . . .
if (nargin < 2)
    bins = 10;
elseif (numel(bins) ~= 1)
    error('Second argument must be a scalar when it is supplied.');
end

% Make equally spaced bin centers
bin_width = 2 * pi / bins;
bin_center = 0 : bin_width : (2*pi - bin_width);

% To make the histogram, we first assign each number to a bin based on its angle.
% We compute (1 - cos(angle(z) - ref_angle)) as a convenient way of comparing the
% angle of each data value with a list of reference angles (this gives us minima
% when the angle difference is a multiple of 2*pi and maxima at odd multiples of pi).
dist_to_bin = (1 - cos(repmat(angle(Z), [1,bins]) - repmat(bin_center, [size(Z,1),1])));

% We find the best bin by subtracting the best fit from each row and then looking for
% zeros; this should label the best fit reference angle for each row with a 1 . . .
bin_scores = (dist_to_bin - repmat(min(dist_to_bin, [], 2), [1, bins]) == 0);

% . . . except for the annoying case where a data point falls between two bins.
% The only fair way to arbitrate this is to choose one randomly; although this
% is not efficient, it avoids the bias of choosing a default (e.g., always go
% closer to 0 degrees).
bin_ambiguity = find(sum(bin_scores, 2) > 1);  % any ambiguous points?
for ambi = 1:length(bin_ambiguity)
    ambiguous = bin_ambiguity(ambi);
    ambins = find(bin_scores(ambiguous,:));   % who are the contenders?
    ambins = ambins(randperm(length(bins)));  % mix 'em up
    bin_scores(ambiguous, :) = 0;             % zero the scores
    bin_scores(ambiguous, ambins(1)) = 1;     % and set our winner
end

% Turn the scores into bin numbers
bin_assign =  bin_scores * (1:bins)';

% Now we add up the magnitudes of all data points in each bin and keep
% this as the weight of the bin.
histogram = full(sparse(bin_assign, 1, abs(Z), bins, 1));

if (nargout > 0)
    weights = histogram';
    if (nargout > 1)
        anglebins = bin_center;
    end
else % we use patches because they look nicer than the lines in 'rose'
    xbars = zeros(3, 0);
    ybars = zeros(3, 0);
    for bin = 1:bins  % draw a triangle for each bin
        [x, y] = pol2cart([0; bin_center(bin)-(bin_width/2); bin_center(bin)+(bin_width/2)], ...
                          [0; histogram(bin); histogram(bin)]);
        xbars = cat(2, xbars, x);
        ybars = cat(2, ybars, y);
    end
    cla;
    polar(0, max(histogram));   % make a backdrop
    patch(xbars, ybars, 'b');   % and draw the patches
end
