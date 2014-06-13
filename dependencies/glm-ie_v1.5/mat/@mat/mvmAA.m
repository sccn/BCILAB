% base matrix: matrix vector multiplication with A.*conj(A)
% Note that (A.*conj(A))*x = diag(A'*diag(x)*A).
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2012 January 05

function y = mvmAA(A,x,ctransp)

  if nargin<3, ctransp = false; end
  if ctransp
    if strcmp(type(A),'vcat')
      y = zeros(size(A,2),1);
      arg = args(A);
      ii = 0;
      for i=1:numel(arg)
        j = size(arg{i},1);
        y = y + vec(mvmAA(arg{i},x(ii+(1:j)),ctransp));
        ii = ii+j;
      end
    else
      fA = full(A)';
      y = (fA.*conj(fA))*x;
    end
  else
    if strcmp(type(A),'vcat')
      y = zeros(size(A,1),1);
      arg = args(A);
      ii = 0;
      for i=1:numel(arg)
        j = size(arg{i},1);
        y(ii+(1:j)) = mvmAA(arg{i},x,ctransp);
        ii = ii+j;
      end  
    else
      fA = full(A);
      y = (fA.*conj(fA))*x;
    end
  end

  
function y = vec(x)
  y = x(:);