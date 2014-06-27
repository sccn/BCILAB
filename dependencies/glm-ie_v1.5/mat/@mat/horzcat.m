% base matrix: horizontal concatenation
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 04

function C = horzcat(varargin)

  if nargin==1, C = varargin{1}; return, end                      % trivial case
  szs = zeros(nargin,2);
  for i=1:nargin
    if ~isa(varargin{i},'mat') && ...
       ~isnumeric(varargin{i}), error('Wrong type'), end
    szs(i,:) = size(varargin{i});
    varargin{i} = [varargin{i}'];                   % transposition of arguments
  end
  if any(diff(szs(:,1))~=0), error('Matrices need same number of rows'), end
  C = mat([sum(szs(:,2)),szs(1,1)],1,varargin{:},'hcat')';       % transposition