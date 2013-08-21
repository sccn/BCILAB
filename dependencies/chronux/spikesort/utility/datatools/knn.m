function [neighbors, distances] = knn(data, k)
%KNN               Finds the K nearest neighbors in a reference matrix.
%   NEIGHBORS = KNN(DATA, K) takes the N x D matrix DATA, where each row
%   is treated as a D-dimensional point, and finds the K nearest neighbors
%   of each point (using a Euclidean distance metric).  NEIGHBORS is an
%   N x K matrix of indices into the rows of the DATA such that the jth
%   nearest neighbor of the ith row of DATA is NEIGHBORS(i,j).
%
%   [NEIGHBORS, DISTS] = KNN(DATA,K) also returns the matrix DISTS of
%   distances such that DISTS(i,j) is the distance:
%         || DATA(i,:) - DATA(NEIGHBORS(i,j),:) ||

% Currently uses a naive O(N^2) algorithm, (although at least
% its broken up to be only O(N) memory).

[N,D] = size(data);
default_chunk = 500;

if (nargin < 2)
	k = 1;
end
% go a chunk at a time so we're not holding the entire N^2 matrix in memory at once
if (N > default_chunk)
    chunksize = default_chunk;
else
    chunksize = N;
end

neighbors = zeros(N, k); %  neighbor indices
distances = zeros(N, k);

for chunkstart = 1:chunksize:N        
    chunkfinis = min(chunkstart+chunksize-1, N);
    disp(['Examining waveforms ' num2str(chunkstart) ' to ' num2str(chunkfinis) '.']);
    
    dists = pairdist(data(chunkstart:chunkfinis,:), data, 'nosqrt', 'reuse');
	[dists, inds] = sort(dists, 2);

    neighbors(chunkstart:chunkfinis, :) = inds(:, 2:(k+1));
	distances(chunkstart:chunkfinis, :) = dists(:, 2:(k+1));
end

if (nargout > 2),  distances = sqrt(distances);  end;
