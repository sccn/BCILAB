% base matrix: generic getter
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 February 14

function field = get(A,fieldname)
  
  A = struct(A); field = [];
  if isfield(A,fieldname)
    eval(['field=A.',fieldname,';']);
  else
    error('Field does not exist.')
  end  
