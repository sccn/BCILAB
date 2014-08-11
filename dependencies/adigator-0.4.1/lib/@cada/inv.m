function y = inv(x)
% CADA overloaded version of INV.
% If x is a matrix, then this function calls cadainversederiv to compute
% the derivatives of the matrix inverse. If x is a scalar, then 1./x is
% called. See cadainversederiv for information on computation of the
% derivative inverse.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end
NUMvod = ADIGATOR.NVAROFDIFF;
PFLAG  = ADIGATOR.PRINT.FLAG;
% --------------------------- Parse Input ------------------------------- %
if x.func.size(1) ~= x.func.size(2)
  error('Can only Invert Square Matrices')
else
  N = x.func.size(1);
end

if N == 1
  y = 1./x;
  return
end

%-------------------------------------------------------------------------%
%                      Build Function Properties                          %
%-------------------------------------------------------------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[N N],'zerolocs',[],...
  'value',[]);

if PFLAG
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
end

if ~isempty(x.func.value)
  y.func.value = inv(x.func.value);
elseif ~isempty(x.func.zerolocs)
  xtemp = rand(N,N);
  xtemp(x.func.zerolocs) = 0;
  warning('off','all')
  ytemp = inv(xtemp);
  warning('on','all')
  y.func.zerolocs = find(~ytemp(:));
  ytemp = abs(ytemp);
else
  ytemp = ones(N,N);
end

% ------------------------Build Derivative Properties----------------------
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    nzlocs = cadainversederiv(x,ytemp,Vcount,derivstr,DPFLAG);
    if ~isempty(nzlocs)
      y.deriv(Vcount).nzlocs = nzlocs;
      y.deriv(Vcount).name   = derivstr;
    end
  end
end

% --------------------------Function Printing --------------------------- %
if PFLAG
  fprintf(fid,[indent,funcstr,' = inv(',x.func.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;