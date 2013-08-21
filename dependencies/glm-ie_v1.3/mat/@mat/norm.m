% base matrix: norm of the full matrix
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 15

function y = norm(A)
  y = norm(full(A));