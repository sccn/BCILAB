% base matrix: convert to full matrix, relies on size and mtimes
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function B = full(A)
  [m,n] = size(A);
  B = zeros(m,n); ei = zeros(n,1);
  for i=1:n, ei(i) = 1; B(:,i) = A*ei; ei(i) = 0; end