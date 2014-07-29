function x = adigatorCreateAuxInput(xsize,varargin)
% Syntax
%  x = adigatorCreateAuxInput(xsize)
% This function creates an input,x, (which has NO derivatives, but is not 
% a fixed value), that is to be used with the function source-to-derivative
% source transformation function, adigator.
%             OR
%  x = adigatorCreateAuxInput(xsize,value)
%  This is useful if you have a vectorized input which is actually a fixed
%  value. If the first dimension is vectorized, then value should be a row
%  vector, if second dimension is vectorized, then value should be a column
%  vector.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also: adigatorCreateDerivInput adigatorOptions adigator

if isequal(size(xsize),[1 2]) && isequal(xsize,floor(xsize))
  func.size = xsize;
else
  error('first input xsize must be integer array of size 1 by 2')
end
if isinf(xsize(1)) && isinf(xsize(2))
  error('only one dimension of the input may be vectorized');
end

if nargin == 2
  value   = varargin{1};
  valsize = size(varargin{1});
  if any(xsize(~isinf(xsize)) ~= valsize(~isinf(xsize))) || any(valsize(isinf(xsize))~=1)
    error('Invalue value input');
  end
  func.value = value;
end

x = cada(0,func,[]);
