% base matrix: Kronecker product
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 29

function C = kron(A,B)

  C = mat(size(A).*size(B),1,A,B,'kron');