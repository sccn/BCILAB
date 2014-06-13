% base matrix: return max element, relies on A.sz
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 12

function e = end(A,id,dummy)

if id>2
  e = 1;
else
  e = size(A,id);
end