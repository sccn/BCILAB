function rcov = cov_blockgeom(X,blocksize)
% like cov(), just robust (using the blockwise geometric median)
% The blocksize allows to reduce the memory requirements, at the cost of reduced robustness
% against outliers that occupy fewer samples than the blocksize.

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

if ~exist('blocksize','var')
    blocksize = 10; end

[n,m] = size(X);
X = bsxfun(@minus,X,median(X));
U = zeros(length(1:blocksize:n),m*m);
for k=1:blocksize
    range = min(n,k:blocksize:(n+k-1));
    U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,m),reshape(X(range,:),[],m,1)),size(U));
end
rcov = real(reshape(geometric_median(U/blocksize),m,m));
