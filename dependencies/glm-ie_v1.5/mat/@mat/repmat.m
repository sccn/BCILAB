% base matrix: replication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 29

function C = repmat(A,r)
  
  r = [1,1].*r(:)';                            % make sure, r has two components
  C = mat(r.*size(A),1,A,r,'rep');