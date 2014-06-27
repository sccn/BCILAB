% base matrix: diagonal of full matrix
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 12

function B = diag(A)

  B = diag(full(A));