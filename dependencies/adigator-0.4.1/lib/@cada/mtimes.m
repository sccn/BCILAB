function z = mtimes(x,y)
% CADA overloaded version of function MTIMES.
% if x and y are matrices, then this function calls cadamtimesderiv to
% compute the derivatives of the matrix multiplication. If x or y is a
% scalar, then z = x.*y is called instead. See cadamtimesderiv for
% information on the computation of matrix inverse.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  z = cadaEmptyEval(x,y);
  return
end
PFLAG  = ADIGATOR.PRINT.FLAG;

NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;

%NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% ----------------------------Parse Inputs------------------------------- %
if isa(x,'cada') && isa(y,'cada')
  % Both Inputs are Symbolic
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  yMrow = y.func.size(1); yNcol = y.func.size(2);
elseif isa(x,'cada')
  % y is numeric input
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  [yMrow,yNcol] = size(y);
  ytemp.id = [];
  ytemp.func = struct('name',[],'size',[yMrow yNcol],'zerolocs',[],...
    'value',y);
  if PFLAG
    if yMrow*yNcol == 1
      ytemp.func.name = num2str(y,16);
    else
      ytemp.func.name = cadamatprint(y);
    end
  end
  ytemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  y = ytemp;
  y = class(y,'cada');
else
  % x is numeric input
  yMrow = y.func.size(1); yNcol = y.func.size(2);
  [xMrow,xNcol] = size(x);
  xtemp.id = [];
  xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
    'value',x);
  if PFLAG
    if xMrow*xNcol == 1
      xtemp.func.name = num2str(x,16);
    else
      xtemp.func.name = cadamatprint(x);
    end
  end
  xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  x = xtemp;
  x = class(x,'cada');
end
% ----------------------------Function Sizing------------------------------
if (xMrow == 1 && xNcol == 1) || (yMrow == 1 && yNcol == 1)
  z = x.*y;
  return
elseif xNcol == yMrow
  FMrow = xMrow;
  FNcol = yNcol;
else
  error('Inputs are not of compatible sizes');
end

%-------------------------------------------------------------------------%
%                      Build Function Properties                          %
%-------------------------------------------------------------------------%
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();

z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);
if isinf(FMrow) && ~isinf(FNcol)
  xMrow = 1; FMrow = 1; zvec = 1;
  if cadaCheckForDerivs(y)
    error(['Cannot perform mtimes with vectorized input times ',...
      'non-vectorized input when the non-vectorized input has derivatives']);
  end
elseif isinf(FNcol) && ~isinf(FMrow)
  yNcol = 1; FNcol = 1; zvec = 2;
  if cadaCheckForDerivs(x)
    error(['Cannot perform mtimes with vectorized input times ',...
      'non-vectorized input when the non-vectorized input has derivatives']);
  end
elseif isinf(FMrow) || isinf(FNcol)
  error(['Result of vectorized mtimes must be of dimension N by n or ',...
    'n by N, where N is vectorized dimension length'])
else
  zvec = 0;
end

% --------------Function Numeric Values and Sparsity--------------------- %
if ~isempty(x.func.value) && ~isempty(y.func.value)
  % z is numeric
  z.func.value = x.func.value*y.func.value;
else
  spflag = 0;
  if ~isempty(x.func.value)
    xtemp  = logical(x.func.value);
    spflag = 1;
  elseif ~isempty(x.func.zerolocs)
    xtemp = true(xMrow,xNcol);
    xtemp(x.func.zerolocs) = false;
    spflag = 1;
  else
    xtemp = true(xMrow,xNcol);
  end
  if ~isempty(y.func.value)
    ytemp  = logical(y.func.value);
    spflag = 1;
  elseif ~isempty(y.func.zerolocs)
    ytemp = true(yMrow,yNcol);
    ytemp(y.func.zerolocs) = false;
    spflag = 1;
  else
    ytemp = true(yMrow,yNcol);
  end
  if spflag == 1
    ztemp = (1*xtemp)*ytemp;
    z.func.zerolocs = find(~ztemp(:));
    if length(z.func.zerolocs) == FMrow*FNcol
      z.func.zerolocs = [];
      z.func.value    = zeros(FMrow,FNcol);
    end
  end
end

% ------------------------Build Derivative Properties----------------------
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod;
  if ~isempty(x.deriv(Vcount).nzlocs) || ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    if zvec
      nzlocs = cadamtimesderivvec(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,zvec);
    else
      nzlocs = cadamtimesderiv(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,'mtimes');
    end
    if ~isempty(nzlocs)
      z.deriv(Vcount).nzlocs = nzlocs;
      z.deriv(Vcount).name   = derivstr;
    end
  end
end
% --------------------------Function Printing --------------------------- %
if PFLAG == 1
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'*',y.func.name,';\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
