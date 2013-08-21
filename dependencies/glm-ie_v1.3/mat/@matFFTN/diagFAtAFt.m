% n-dimensional FFT matrix: (real) diagonal in Fourier space of A'*A
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 20

function d = diagFAtAFt(A)

  d = ones(A.imsz);