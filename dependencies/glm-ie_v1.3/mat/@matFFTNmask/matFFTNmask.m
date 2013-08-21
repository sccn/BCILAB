% n-dimensional partial FFT matrix: constructor
%
% A = matFFTNmask(mask,complex,full)
%
% The transform is centered.
%
%  mask    [m,n,..]: indices to keep
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%  full    [1,1]: flag saying whether A*x contains removed indices as zeros or
%                 only the non-removed parts,                      default false
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 23

function A = matFFTNmask(mask,complex,full)

  A.mask = logical(mask); A.imsz = size(A.mask);                    % image size
  sz = [sum(A.mask(:)),prod(A.imsz)];                              % matrix size
  if nargin<2,  complex = [1,1]; end                               % set default
  if nargin>2,  A.full = full; else A.full = false; end

  if A.full, sz(1) = sz(2); end
  A = class(A,mfilename,mat(sz,complex));