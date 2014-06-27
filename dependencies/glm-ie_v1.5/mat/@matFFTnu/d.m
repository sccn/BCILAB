% Non-uniformly spaced 2d/3d FFT matrix: derivative
%
% (c) by Hannes Nickisch & Alexander Loktyushin, 
%                               MPI for Biological Cybernetics, 2011 February 12

function dAi = d(A,i)

  if nargin<2                                         % default derivative index
    i = 0;
  else
    dimax = A.ndims;                                  % maximal derivative index
    i = fix(i);
    if i>dimax, error('i is too large.'), end
    if i<0,     error('i is too small.'), end

    if strcmp(A.type,'fessler')
      if numel(A.dP{i})==0                                          % precompute
        [junk, A.dP{i}] = nufft_init(A.k, A.imsz, A.fftsz, A.kbwsz, i);
      end
    elseif strcmp(A.type,'ferrara')
      error('Derivative for ferrara not implemented.')
    end
  end
  A.di = i;                            % which derivative, 0 means no derivative
  dAi = A;                                                       % assign output