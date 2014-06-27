% Non-uniformly spaced 2d/3d FFT matrix: matrix vector multiplication with A'*A
%                                        which is symmetric block Toeplitz
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                                  MPI for Biological Cybernetics, 2011 March 03

function y = mvmAtA(A,x)

if any([strfind(lower(A.type),'near'), strfind(lower(A.type),'lin'), ...
                                             strfind(lower(A.type),'cub')])
  k = bsxfun(@minus, A.k, A.imsz(:)+1 );                       % subtract center
  k = bsxfun(@times, k, 1./(2*A.imsz(:)) );          % get back scale, padding 2
  c = A.f*2*ones(size(k,2),1);                                 % constant vector
  w = A.P'*c;         
  % w = matFFTN(2*A.imsz)*(matFFTnu(2*A.imsz,k',A.type)'*c);
  y = A.F'*(w.*(A.F*x));
else
  y = []; error('Only implemented for types nearest, linear and cubic.');
end