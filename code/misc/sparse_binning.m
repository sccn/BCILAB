function [B,R] = sparse_binning(V,rows,columns)
% Round and bin the given values into a sparse array.
% [BinIndices,Residuals] = sparse_binning(Values,Padding)
% 
% In:
%   Values : vector of values that yield valid indices (>1) when rounded.
%
%   Rows : number of rows to reserve (default: []=as many as needed)
%
%   Columns : number of columns to reserve (default: []=as many as needed)
%
% Out:
%   BinIndices : sparse array containing indices into Values, where the horizontal axis is the
%                rounded value and the vertical axis are the ranks of values that fall into the same
%                bin.
%
%   Residuals : optionally the fractional offset between the values and the bin indices.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-12

[V,order] = sort(V(:)'); %#ok<TRSRT>
% get the bin index where each value falls
bins = round(V);
% get the within-bin rank of each value
ranks = [false ~diff(bins)];
ranks = 1+cumsum(ranks)-cummax(~ranks.*cumsum(ranks));

% create a sparse matrix whose k'th column holds the indices to values in the k'th bin
if nargin == 1
    B = sparse(ranks,bins,order);
elseif nargin == 2
    B = sparse(ranks,bins,order,rows,max(bins));
else
    if isempty(rows)
        rows = max(ranks); end
    if isempty(columns)
        columns = max(bins); end
    B = sparse(ranks,bins,order,rows,columns);
end

if nargout>0
    R = V-bins; end