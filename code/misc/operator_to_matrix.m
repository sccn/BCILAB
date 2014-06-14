function M = operator_to_matrix(op,n)
% Calculates a transformation matrix from a linear operator (expressed as a function).
% Matrix = operator_to_matrix(Operator,NumDimensions)
%
% This function allows to specify linear operators using concise MATLAB syntax (e.g., @(x)diff(x'))
% and can generate the necessary matrix from it that is used in numerical solvers, etc.
%
% The output is cached by default.
%
% In:
%   Operator : the linear operator to apply.
%
%   NumDimensions : the number of dimensions of the input vectors that the operator applies to.
%
% Out:
%   Matrix : a sparse matrix that captures the action of the operator.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-02-04

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

persistent operator_cache;
opid = ['x' hlp_cryptohash({op,n})];
try
    M = operator_cache.(opid);
catch
    % convert a linear operator function to a matrix (given the input dimensionality)
    % this works by going through the canonical basis vectors of the space and projecting them
    % one-by-one through the operator
    fprintf('Evaluating and caching operator...');
    M = hlp_diskcache('general',@operator_to_matrix_cached,op,n);
    operator_cache.(opid) = M;
    fprintf('done.\n');
end

function M = operator_to_matrix_cached(op,n)
vec = @(x)x(:);
M = sparse([]);
w = zeros(n,1);
for c=[n 1:n-1]
    v = w; v(c) = 1;
    M(:,c) = vec(op(v));
end
