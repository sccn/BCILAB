% 2d wavelet transform matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 19

function y = mvm(A,x,ctransp)
  x = reshape(x,A.imsz);
  if ctransp
    if isreal(x)
      y = fwtn(x,A.nl,A.qmf,1);
    else
      y = fwtn(real(x),A.nl,A.qmf,1) + 1i*fwtn(imag(x),A.nl,A.qmf,1);
    end
  else
    if isreal(x)
      y  = fwtn(x, A.nl, A.qmf);
    else
      y  = fwtn(real(x), A.nl, A.qmf) + 1i*fwtn(imag(x), A.nl, A.qmf);
    end
  end