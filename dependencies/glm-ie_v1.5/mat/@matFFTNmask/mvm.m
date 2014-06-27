% n-dimensional partial FFT matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 20

function y = mvm(A,x,ctransp)
  
  N = prod(A.imsz);
  if ctransp
    if A.full
      y = reshape(x,A.imsz);
    else
      y = zeros(A.imsz);                            % allocate memory for output
      y(A.mask) = x;
    end
    if numel(A.imsz)<3
      y = ifft2(ifftshift(y))*sqrt(N);
    else
      y = ifftn(ifftshift(y))*sqrt(N);
    end
  else
    if numel(A.imsz)<3
      y = fftshift( fft2( reshape(x,A.imsz) )/sqrt(N) );
    else
      y = fftshift( fftn( reshape(x,A.imsz) )/sqrt(N) );  
    end
    if A.full
      y(~A.mask) = 0;
    else
      y = y(A.mask);
    end
  end