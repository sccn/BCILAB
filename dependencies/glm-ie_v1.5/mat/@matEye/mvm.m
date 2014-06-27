% Identity matrix example
function y = mvm(A,x,ctransp)
  sz = size(A);                                                           % size
  p = A.par;                                              % additional parameter
  if ctransp                                                        % transposed
    fprintf('We did an MVM with eye(%d,%d)'' par=%s.\n',sz,p)
    y = x;
  else                                                          % not transposed
    fprintf('We did an MVM with eye(%d,%d) par=%s.\n',sz,p)
    y = x;
  end