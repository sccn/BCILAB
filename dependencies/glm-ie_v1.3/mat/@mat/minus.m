% base matrix: minus
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 23

function C = minus(A,B)

  sa = size(A); sb = size(B); if any(sa~=sb), error('Wrong matrix sizes'), end
  C = mat(sa,1,A,-B,'sum');