function markers = leadingedges(data)
%LEADINGEDGES      Marks 0 -> NONZERO transitions along data columns.
%   EDGEMARKERS = LEADINGEDGES(DATA) takes an M x N matrix DATA and
%   returns a M x N matrix EDGEMARKERS, containing a '1' in each location
%   that corresponds to the first nonzero value in a DATA column following
%   one or more zeros above it in the same column.  If the first row of
%   any DATA column contains a nonzero value, it will be marked in
%   EDGEMARKERS.
%
%   DATA can be a matrix of type double, uint8, or logical; in all of
%   these cases, EDGEMARKERS is of type logical.  DATA should not contain
%   NaN values.

if (issparse(data)),  error('Sparse matrices can not be used with LEADINGEDGES.');  end

markers = CORE_leadingedges(data);
