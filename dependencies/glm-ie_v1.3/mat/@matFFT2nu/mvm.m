% Non-uniformly spaced 2d FFT matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 17

function y = mvm(A,x,ctransp)
  
  N = prod(A.imsz);
  if ctransp
    Ptx = reshape(full(A.P'*x), A.fftsz);  % apply sparse interpolation matrix P
    y = prod(A.fftsz)/sqrt(N) * ifft2(Ptx);              % apply oversampled fft
    y = y(1:A.imsz(1),1:A.imsz(2),:);                   % eliminate zero padding
    y = conj(A.d) .* y;                        % apply diagonal scaling matrix D
  else
    Sx = A.d .* reshape(x,A.imsz);             % apply diagonal scaling matrix D
    FSx = fft2(Sx,A.fftsz(1),A.fftsz(2))/sqrt(N);        % apply oversampled fft
    y = A.P*FSx(:);                        % apply sparse interpolation matrix P
  end