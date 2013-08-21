function dim = isvectord(vect)
%ISVECTORD         Returns the orientation of a 1-D vector.
%   ISVECTORD(VECT) is non-zero if exactly one dimension of VECT has
%   length greater than 1.  The return value is then the index of that
%   dimension.  Note that NDIMS can not be used to decide this question,
%   because it returns 2 for, e.g., (M x 1) and (1 x M) arrays.
%
%   Example:
%      isvectord(1);             % returns 0
%      isvectord([1 2 ; 3 4])    % returns 0
%      isvectord([1:10])         % returns 2
%
%   See also ISVECTOR (Matlab R14 and later).

nonsingle = (size(vect) > 1);
dim = find(nonsingle);
if ((length(dim)>1) || (isempty(dim))),  dim = 0;  end;