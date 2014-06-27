% base matrix: matrix vector multiplication with A.*conj(A)
% Note that (A.*conj(A))*x = diag(A'*diag(x)*A).
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2012 August 18

function y = mvmAA(A,x,ctransp)
  x = reshape(x,A.imsz); if nargin<3, ctransp=0; end
  if length(A.qmf)~=2 || norm(A.qmf(:)-[1;1]/sqrt(2))>1e-12
    error('Implementation for Haar wavelets only.')
  end

  if ctransp
    if isreal(x)
      y = fwtn(x,A.nl,A.qmf,1,1);
    else
      y = fwtn(real(x),A.nl,A.qmf,1,1) + 1i*fwtn(imag(x),A.nl,A.qmf,1,1);
    end
  else
    if isreal(x)
      y  = fwtn(x,A.nl,A.qmf,0,1);
    else
      y  = fwtn(real(x),A.nl,A.qmf,0,1) + 1i*fwtn(imag(x),A.nl,A.qmf,0,1);
    end
  end