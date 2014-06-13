% 2d convolution matrix: (real) diagonal in Fourier space of A'*A
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function d = diagFAtAFt(A,opt)

  if strcmp(A.shape,'circ') || (nargin>1 && strcmp(opt,'cheap'))
    d = abs(A.F); d = fftshift(d.*d);
  else
    sz = size(A.F);
    F  = matFFTNmask(ones(sz));
    d = real(reshape(diag(F*A'*A*F'),sz));
  end
