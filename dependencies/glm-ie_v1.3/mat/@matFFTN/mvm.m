% n-dimensional FFT matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 January 31

function y = mvm(A,x,ctransp)

  y = reshape(x,A.imsz);
  N = prod(A.imsz(A.dims));
  if ctransp
    for d=A.dims, y=fftshift(ifft(ifftshift(y,d),[],d),d); end, y = y*sqrt(N);
  else
    for d=A.dims, y=fftshift( fft(ifftshift(y,d),[],d),d); end, y = y/sqrt(N);
  end