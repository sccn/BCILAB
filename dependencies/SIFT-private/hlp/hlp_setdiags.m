function C = hlp_setdiags(C,vals)
% Assign values to the diagonals of a multidimensional matrix
% Input: 
%   C:      A multidimensional matrix. Dimensions of C greater than two are
%           considered pages of 2D matrices whose diagonals will be modified.
%   vals:   Value(s) to store on the diagonal. This can be a scalar or a 
%           vector containing values to sequentially assign to diagonals in 
%           column-major order.
% 
% Example 1: zero out the diagonals of a 4D matrix
% C = rand(3,3,5,6);
% C = hlp_setdiags(C,0)
% 
% Example 2: set the diagonals of a 3D matrix to an ascending int sequence:
% C = rand(3,3,5);
% C = hlp_setdiags(C, [1:3*5])
% 
% Author: Tim Mullen, SCCN/INC/UCSD, 10-5-2013

[d1, d2, d3] = size(C); 

if ~isscalar(vals) && numel(vals)~=d1*d3
    error('vals must either be scalar or its length must equal the total number of diagonal elements in C');
end

ind = vec(bsxfun(@plus,(1:d1+1:d1*d2)',d1*d2*(0:d3-1)));
C(ind) = vals;