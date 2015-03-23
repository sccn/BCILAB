function indices = utl_equidistant_subset(x,y,z,num,projected)
% Find a subset of 3d points that are maximally evenly distributed.
% Indices = utl_equidistant_subset(X,Y,Z,SubsetSize,Projected)
%
% This uses a greedy algorithm which starts with the empty set and then successively adds the
% maximally distant point to the set.
%
% In:
%   X : vector of X coordinates
%
%   Y : vector of Y coordinates
%
%   Z : vector of Z coordinates
%
%   SubsetSize : number of elements in the desired subset
%
%   Projected : whether the points shall be projected onto a sphere centered on the coordinate
%               origin (default: true)
%
% Out:
%   Indices : indices of the selected points
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-06-22

if ~exist('projected','var') || isempty(projected)
    projected = true; end

x = x(:);
y = y(:);
z = z(:);

% selected indices so far
indices = 1;
while length(indices) < num
    % candidate indices
    candidates = setdiff(1:length(x),indices);
    % positions in current set
    in_set = [x(indices) y(indices) z(indices)];
    % positions not in current set
    out_set = [x(candidates) y(candidates) z(candidates)];
    % project positions onto sphere
    if projected
        in_set = bsxfun(@times,in_set,1./sqrt(sum(in_set.^2,2)));
        out_set = bsxfun(@times,out_set,1./sqrt(sum(out_set.^2,2)));
    end
    % calc pairwise distances
    distances = sqrt(sum(bsxfun(@minus,permute(in_set,[3 1 2]),permute(out_set,[1 3 2])).^2,3));
    % get minimum of that to get distance of each candidate to the in-set
    set_distances = min(distances,[],2);
    % find candidate with max distance to in-set
    [dummy,cand_index] = max(set_distances); %#ok<ASGLU>
    indices = [indices candidates(cand_index)]; %#ok<AGROW>
end

indices = sort(indices); 
