% 2d partial FFT matrix (orthonormal) selecting rows: constructor
%
% A = matFFT2line(imsz,id,complex)
%
%  sz      [1,2]: size of the image to transform
%  id      [1,p]: row indices to keep, default id=1:sz(1); the transform is
%                 centered i.e. floor(sz(1)/2+1) is the center line
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%  full    [1,1]: flag saying whether A*x contains removed indices as zeros or
%                 only the non-removed parts,                      default false
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function A = matFFT2line(imsz,id,complex,full)

  A.imsz = imsz; if numel(A.imsz) == 1, A.imsz = A.imsz*[1,1]; end  % image size
  if nargin>1
    id(id<1) = 1; id(id>A.imsz(1)) = A.imsz(1);              % put into interval
    A.id = unique(fix(id(:)));              % subset of 1..imsz(1), double is ok
  else
    A.id = 1:A.imsz(1);
  end
  sz = [numel(A.id)*A.imsz(2),prod(A.imsz)];                       % matrix size
  if nargin<3,  complex = [1,1]; end                               % set default
  if nargin>3,  A.full = full; else A.full = false; end
  if A.full, sz(1) = sz(2); end
  A = class(A,mfilename,mat(sz,complex));
