% wavelet transform matrix: constructor
%
% A = matWav(sz,qmf,nl,complex)
%
%  sz      [1,2]: size of the tensor to transform
%  qmf     [1,p]: quadrature mirror filter, default is Haar qmf = [1,1]/sqrt(2)
%  nl      [1,1]: number of layers,         default is full nl
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 26

function A = matWav(sz,qmf,nl,complex)

  imsz = sz(:)';  A.imsz = imsz;                                    % image size
  if numel(A.imsz)==1, A.imsz = [A.imsz,1]; end
  sz = prod(imsz)*[1,1];                                           % matrix size
  if nargin<4, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  % quadrature mirror filter
  A.qmf = [1,1]/sqrt(2);                                       % default is Haar
  if nargin>1, if numel(qmf)>0, A.qmf = qmf; end, end
  imsz = imsz/2;                              % compute maximum number of layers
  for nlmax = 1:max(min(fix(log2(imsz))),1)
    if any(imsz ~= fix(imsz)), nlmax = nlmax-1; break, end
    imsz = imsz/2;
  end
  A.nl = nlmax;                                    % default is maximum possible  
  if nargin>2, if numel(nl)>0, A.nl = min(fix(nl),nlmax); end, end  
  A = class(A,mfilename,mat(sz,complex));