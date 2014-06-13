% 2d convolution matrix: (real) diagonal in Fourier space of A'*A
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 12

function d = diagFAtAFt(A)

  if strcmp(A.shape,'circ')
    d = abs(A.F); d = fftshift(d.*d);
  else
    sz = size(A.F);
    F  = matFFTNmask(ones(sz));
    d = real(reshape(diag(F*A'*A*F'),sz));
  end