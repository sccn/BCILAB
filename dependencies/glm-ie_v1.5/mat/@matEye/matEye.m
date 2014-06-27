% Identity matrix example
function A = matEye(n,m, par, complex)
  sz = [n,m];                                               % size of the matrix
  if nargin<3, A.par = ''; else A.par = par; end         % additional parameters
  if nargin<4, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));