% base matrix: display something, relies on size
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 22

function disp(A)

  str = class(A);
  if strcmp(str,'mat'), str = [str,'_',A.type]; end
  fprintf('<%dx%d %s>\n', size(A,1), size(A,2), str);