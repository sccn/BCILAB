function varargout = adigatorProjectVectLocs(N,varargin)
% Syntax for adigatorProjectVectLocs is
% [I1,I2,...,In] = adigatorProjectVectLocs(N,i1,i2,...,in)
% Where N is vectorized dimension,
% i1,...,in are the locations of the small problem, and
% I1,...,In are the locations of the large problem
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
if nargout ~= nargin-1 || nargin == 1 
error('invalid inputs')
end
varargout = cell(1,nargout);
for ii = 1:nargout
  n = length(varargin{ii});
  varargout{ii} = repmat((1:N).',1,n) + N.*ones(N,n)*diag(varargin{ii}-1);
end
