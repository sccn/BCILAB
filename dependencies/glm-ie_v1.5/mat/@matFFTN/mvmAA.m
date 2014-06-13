% base matrix: matrix vector multiplication with A.*conj(A)
% Note that (A.*conj(A))*x = diag(A'*diag(x)*A).
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2012 January 05

function y = mvmAA(A,x,ctransp)

  N = prod(A.imsz(A.dims));
  y = ones(A.imsz)*sum(x(:))/N;