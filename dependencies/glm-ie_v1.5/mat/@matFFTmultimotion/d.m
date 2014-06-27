%    derivative

function dAi = d(A,i)

  if nargin<2                                         % default derivative index
    i = 0;
  else
    dimax = A.num_patches;                       % maximal derivative index
    i = fix(i);
    if i>dimax, error('i is too large.'), end
    if i<1,     error('i is too small.'), end
  end                             
  A.di = i;               % which derivative, 0 means all of the respective kind
  dAi = A;                                                  % assign output
  
  