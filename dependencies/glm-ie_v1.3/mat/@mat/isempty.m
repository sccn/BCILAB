% base matrix: say whether matrix is empty, relies on size
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function y = isempty(A)
  [m,n] = size(A);
  y = m*n==0;