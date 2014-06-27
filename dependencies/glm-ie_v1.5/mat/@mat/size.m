% base matrix: return size, relies on A.sz
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 18

function [sz1,sz2] = size(A,id)
  sz = A.sz;
  if A.ctransp, sz = fliplr(sz); end
  if nargin==1
    if nargout==2
      sz1 = sz(1); sz2 = sz(2);
    else
      sz1 = sz;
    end
  elseif nargin==2
    if nargout>1, error('too many output arguments'), end
    if id(1)>numel(sz)
      sz1 = 1;
    else
      sz1 = sz(id);
    end
  else
    error('wrong number of input arguments')
  end