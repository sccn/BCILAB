% 2d finite difference matrix: constructor
%
% A = matFD2(sz,shape,complex)
%
%  sz      [1,2]: size of the image to transform
%  shape   one of 'zero', 'circ', 'same' or 'valid'               default 'circ'
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function  A = matFD2(sz, shape, complex)

  if any(sz<2), error('[%s.m] size missmatch', mfilename), end
  imsz = sz(:)';  A.imsz = imsz;                                    % image size
  if numel(A.imsz)==1, A.imsz = [A.imsz,1]; end
  if nargin<2, shape = 'circ'; end, A.shape = shape;
  switch shape
    case 'zero'
      sz = prod(imsz)*[2,1];
    case 'circ'
      sz = prod(imsz)*[2,1];
    case 'valid'
      sz = [prod(imsz-[1,0])+prod(imsz-[0,1]), prod(imsz)];
    case 'same'
      sz = prod(imsz)*[2,1];
    otherwise
      error('[matFD2.m] Unknown shape. Choose between ''circ'', ''same'' and ''valid''');
  end

  % A.d is used in diagFAtAFt
  f1 = fftn([1;-1],imsz);
  f2 = fftn([1,-1],imsz);
  A.d = real(fftshift(f1.*conj(f1) + f2.*conj(f2)));

  if nargin<3, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));
