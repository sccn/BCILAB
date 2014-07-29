function y = cadaRemoveRowsCols(x,ySize)
%function y = cadaRemoveRowsCols(x,ySize)
% This function is written to remove any rows/columns from an overloaded
% object. This is called from overloaded binary operations when in the
% Printing run of a FOR loop and two inputs to the binary operation are
% overloaded.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
NUMvod = ADIGATOR.NVAROFDIFF;
NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
y = x;
y.func.size = ySize;
funcstr = ['cada',NDstr,'tempf1'];
y.func.name = funcstr;
ySize(isinf(ySize)) = 1;
x.func.size(isinf(x.func.size)) = 1;
% Function Numerics/Sparsity
if ~isempty(y.func.value)
  y.func.value = y.func.value(1:ySize(1),1:ySize(2));
elseif ~isempty(y.func.zerolocs)
  xtemp = true(x.func.size);
  xtemp(y.func.zerolocs) = false;
  ytemp = xtemp(1:ySize(1),1:ySize(2));
  y.func.zerolocs = find(~ytemp(:));
end
% Get Linear Index
xtemp    = zeros(x.func.size);
xtemp(:) = 1:x.func.size(1)*x.func.size(2);
yref     = xtemp(1:ySize(1),1:ySize(2));
yref     = yref(:);


% Derivatives
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    xrows = x.deriv(Vcount).nzlocs(:,1);
    xcols = x.deriv(Vcount).nzlocs(:,2);
    nzx   = length(xrows);
    dx = sparse(xrows,xcols,1:nzx,prod(x.func.size),ADIGATOR.VAROFDIFF(Vcount).usize);
    dy = dx(yref,:);
    
    [yrows,ycols,derInds] = find(dy);
    if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
    if ~isempty(yrows)
      TFind1 = cadaindprint(derInds(:));
      fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',TFind1,');\n']);
      y.deriv(Vcount).nzlocs = [yrows ycols];
      y.deriv(Vcount).name   = derivstr;
    else
      y.deriv(Vcount).nzlocs = [];
      y.deriv(Vcount).name   = [];
    end
  end
end

% Print Function
if ySize < x.func.size
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'(1:%1.0f,1:%1.0f);\n'],ySize(1),ySize(2));
elseif ySize(1) < x.func.size(1)
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'(1:%1.0f,:);\n'],ySize(1));
else
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'(:,1:%1.0f);\n'],ySize(1));
end