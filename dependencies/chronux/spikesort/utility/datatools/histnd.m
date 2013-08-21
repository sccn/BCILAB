function [ndcounts, inds] = histnd(data, bins)
%HISTND            N-Dimensional histogram
%   [NDCOUNTS, INDS] = HISTND(DATA, BINS);
%   The DATA matrix is treated as a collection of D-dimensional row
%   vectors and the entire range of the data is divided into BINS evenly
%   spaced segments (default: 10) in and the density evaluated on a BINS^D
%   grid.  Rows with any NaN values are ignored.
%
%   The counts of points on this grid are returned as a D-dimensional
%   matrix.  The optional second output containing the bin centers along
%   any axis.

if (nargin == 1)   % default values
	bins = 10;
end

% Make the histogram grid
D = size(data, 2);
ndcounts = zeros(repmat(bins, 1, D));

% We rescale and round the data to break it into (integer numbered) bins
[data, oldmin, oldmax] = rescale(data, 0.5, bins+0.5);
data = round(data);
data(data == 0) = 1;
data(data > bins) = bins;
data = data(~any(isnan(data), 2),:);

% Now compute the density by treating each row of the scaled data as a coordinate in
% the counts matrix and using sparse to do the histogramming.
data = num2cell(data, 1);                 % convert coordinates to cells ...
data = sub2ind(size(ndcounts), data{:});  % and then to unique tags.
ndcounts = reshape(full(sparse(data, 1, 1, bins^D, 1)), size(ndcounts));

if (nargout > 1)
	increment = (oldmax-oldmin)/bins;   % bin width is uniform division of data range
	inds = (increment * (1:bins)) - increment/2 + oldmin;  % compute bin centers
end
