% diagonal matrix (analogue to Matlab's diag)
%
% A = matDiag(a)
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 23

function A = matDiag(a, complex)
  
  sz = numel(a)*[1,1]; A.a = a(:);
  if nargin<2, complex = 1; end
  if numel(complex)==1, complex = complex*[1,1]; end
  A = class(A,mfilename,mat(sz,complex));
