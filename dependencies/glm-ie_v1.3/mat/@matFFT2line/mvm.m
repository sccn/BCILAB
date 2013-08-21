% 2d partial FFT matrix selecting rows: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function y = mvm(A,x,ctransp)
  
  N = prod(A.imsz);
  if ctransp
    if A.full
      y = reshape(x,A.imsz);
    else
      y = zeros(A.imsz);                            % allocate memory for output
      y(A.id,:) = reshape(x,[numel(A.id),A.imsz(2)]);
    end
    y = ifft2(ifftshift(y))*sqrt(N);
  else
    y = fftshift( fft2( reshape(x,A.imsz) )/sqrt(N) );   
    if A.full
      zid = true(A.imsz(1),1); zid(A.id) = 0;
      y(zid,:) = 0;
    else
      y = y(A.id,:);
    end
  end