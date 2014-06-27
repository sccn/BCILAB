% base matrix: type of the matrix
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 12

function str = type(A)

  str = A.type;
  if strcmp(str,'zero')
    str = class(A);
    if strcmp(str,'mat'), str = 'zero'; end
  end