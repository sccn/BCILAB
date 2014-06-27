% base matrix: complex conjugate transpose
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function A = ctranspose(A)
  A.ctransp = xor(A.ctransp,1);