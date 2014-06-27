% base matrix: matrix vector multiplication with A.*conj(A)
% Note that (A.*conj(A))*x = diag(A'*diag(x)*A).
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2012 January 06

function y = mvmAA(A,x,ctransp)
  
  if strcmp(A.type,'conv')
    AA = matConv2(A.f.*conj(A.f), A.sx, A.shape, A.type, A.complex);

    if nargin<3, ctransp = false; end
    if ctransp
      y = AA'*x;
    else
      y = AA*x;
    end
  else
    error('[%s.m] works for ''conv'' only.', mfilename)
  end