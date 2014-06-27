% n-dimensional FFT matrix: constructor
%
% A = matFFTN(sz,dims,complex)
%
% The transform is centered.
%
%  sz      [1,D]: size the data to process
%  dims    [1,d]: dimensions along which the FFT is computed       default 1:D
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 21

function A = matFFTN(sz,dims,complex)

  if numel(sz)<2, sz = [sz(:); 1]; end, A.imsz = fix(sz(:)');       % image size
  sz = prod(sz)*[1,1];                                             % matrix size
  
  if nargin<2, dims = 1:numel(A.imsz); end
  A.dims = unique(fix(dims(:)))';
  if nargin<3,  complex = [1,1]; end                               % set default

  A = class(A,mfilename,mat(sz,complex));