% 2d finite difference matrix: (real) diagonal in Fourier space of A'*A
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 11

function d = diagFAtAFt(A)

  sz = A.imsz;
  if strcmp(A.shape,'circ')
    d = A.d;
  else
    F  = matFFTNmask(ones(sz));
    d  = real(reshape(diag(F*A'*A*F'),sz));
  end
