function m = hlp_minor(X, row, col)
%
% calculate the ijth minor of matrix X by returning the determinant of X
% after removing the ith row and jth column
%
% Input:  
%
%   X:      2- or 3-dimensional matrix
%   row,col to remove
%
% output:  
%
%   m:      if X is 2-D, m = minor(X)
%           if X is 3-D, m(i) = minor(X(:,:,i))
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% 
% 
% Author: Tim Mullen Dec 1st, 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


nd = ndims(X);
if nd < 2 || nd > 3
    error('X must have 2 or 3 dimensions');
end

X(row,:,:) = [];
X(:,col,:) = [];

if nd==2
    m = det(X);
else
    m=zeros(size(X,3),1);
    for k=1:size(X,3)
        m(k) = det(X(:,:,k));
    end 
end

    
