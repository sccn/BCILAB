% base matrix: number of elements in matrix: relies on size
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function N = numel(A,varargin)

  if nargin>1, N=1; return, end 
  [m,n] = size(A); N = n*m;