% base matrix: (real) diagonal in Fourier space of A'*A
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 12

function d = diagFAtAFt(A)

  d = [];
  if strcmp(type(A),'vcat')
    d = 0;
    arg = args(A);
    for i=1:numel(arg)
      di = diagFAtAFt(arg{i});
      if numel(di)>0
        d = d + di;
      else
        d = []; return
      end
    end
  elseif strcmp(type(A),'scale')
    arg = args(A);
    d = arg{1}*arg{1}*diagFAtAFt(arg{2});
  else
    d = [];
  end