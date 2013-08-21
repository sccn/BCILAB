% diagonal matrix: matrix vector multiplication
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 20

function y = mvm(A,x,ctransp)

  if ctransp
    y = conj(A.a(:)) .* x(:);
  else
    y =      A.a(:)  .* x(:);
  end