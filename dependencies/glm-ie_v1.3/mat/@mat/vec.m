% base matrix: vectorisation
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 29

function v = vec(A)

  id.type = '()'; id.subs{1} = ':';
  v = subsref(A,id);