% base matrix: uminus
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 04

function C = uminus(A)

  if size(A,1)==1, C = A*(-1); else C = (-1)*A; end