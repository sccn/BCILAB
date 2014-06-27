% base matrix: vertical concatenation
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 04

function C = vertcat(varargin)

  if nargin<2, error('We need at least 2 inputs'), end
  szs = zeros(nargin,2);
  for i=1:nargin
    if ~isa(varargin{i},'mat') && ...
       ~isnumeric(varargin{i}), error('Wrong type'), end
    szs(i,:) = size(varargin{i});
  end  
  if any(diff(szs(:,2))~=0), error('Matrices need same number of columns'), end
  C = mat([sum(szs(:,1)),szs(1,2)],1,varargin{:},'vcat');