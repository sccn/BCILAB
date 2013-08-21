% base matrix: constructor
%
%  sz      [1,2]: size of the matrix, sz = [m,n]                   default [1,1]
%  complex [1,2]: complex number treatment A:IN^n->OUT^m,          default [1,1]
%                 first number specifies output field OUT^m and the second the
%                 input field IN^m of the matrix A
%                 1 -  real(z) + 1i*imag(z)    IC:   native complex numbers
%                 2 - [real(z); imag(z)]       IR^2: stacked real numbers
%                 3    real(z)                 IR:   real part only
%  A1      [m,n]: matrix of size sz,                         default  zeros(m,n)
%  type   string: {'','hcat','prod','scale','sub','sum','vcat'},      default ''
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 04

function A = mat(sz,complex,varargin)

  if nargin<3
    type = 'zero';         args = [];
    if nargin==0, sz = [1,1]; else if numel(sz)==0, sz = [1,1];  end, end
  elseif nargin==3
    type = 'full';         args = varargin{1};       sz = size(args);
  else
    type = varargin{end};  args = varargin(1:end-1);
  end
  if nargin<2, complex = 1; else if numel(complex)==0, complex = 1; end, end
  if numel(complex)==1, complex = complex*[1,1]; end
  if complex(1)==2, sz(1) = 2*sz(1); end
  if complex(2)==2, sz(2) = 2*sz(2); end
  A.args = args; A.type = type; A.sz = sz; A.complex = complex; A.ctransp = 0;
  A = class(A,'mat');