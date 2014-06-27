function [index_a,index_b] = utl_match_channels(locs_a,locs_b)
% Finds pairs of channels in two given chanlocs structures that are closest.
% [SubsetA,SubsetB] = utl_match_channels(LocationsA,LocationsB)
%
% In:
%   LocationsA : Chanlocs structure or cell array of labels
%
%   LocationsB : Chanlocs structure or cell array of labels
%
% Out:
%   IndexA : a vector of ascending indices into LocationsA; if LocationsA
%            was longer than LocationsB, this will be a subset of channels
%
%   IndexB : a vector of indices into LocationsB ordered so as to pair
%            with the closest matches to LocationsA
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17


% map from labels to full chanlocs
if iscell(locs_a)
    locs_a = hlp_microcache('matchchan',@set_infer_chanlocs,locs_a); end
if iscell(locs_b)
    locs_b = hlp_microcache('matchchan',@set_infer_chanlocs,locs_b); end

% sanity checks
if ~all(isfield(locs_a,{'X','Y','Z'}))
    error('LocationsA does not have associated X,Y,Z coordinate fields.'); end
if ~all(isfield(locs_b,{'X','Y','Z'}))
    error('LocationsB does not have associated X,Y,Z coordinate fields.'); end
if length(unique({locs_a.labels})) < length(locs_a)
    error('LocationsA has redundant labels.'); end
if length(unique({locs_b.labels})) < length(locs_b)
    error('LocationsB has redundant labels.'); end

% strip all locations that have no positions
locs_a(cellfun('isempty',{locs_a.X}) | cellfun('isempty',{locs_a.Y}) | cellfun('isempty',{locs_a.Z})) = [];
locs_b(cellfun('isempty',{locs_b.X}) | cellfun('isempty',{locs_b.Y}) | cellfun('isempty',{locs_b.Z})) = [];

% find closest matches
pos_a = [[locs_a.X]' [locs_a.Y]' [locs_a.Z]'];
pos_b = [[locs_b.X]' [locs_b.Y]' [locs_b.Z]'];
distances = sqrt(sum(bsxfun(@minus,permute(pos_b,[3 1 2]),permute(pos_a,[1 3 2])).^2,3));
[sorted,matches] = sort(distances(:),'ascend'); %#ok<ASGLU>

% generate matched label lists
map = reshape(1:numel(distances),size(distances));
[index_a,index_b] = deal(zeros(1,min(size(distances))));
n=1;
for k=1:length(index_a)
    while ~any(map == matches(n))
        n = n+1; end
    [a,b] = find(map == matches(n),1);
    map(a,:) = 0; map(:,b) = 0;
    index_a(k) = a;
    index_b(k) = b;
    n = n+1;
end

[index_a,ordering] = sort(index_a,'ascend');
index_b = index_b(ordering);
