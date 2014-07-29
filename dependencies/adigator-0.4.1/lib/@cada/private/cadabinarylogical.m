function z = cadabinarylogical(x,y,callerstr)
% Called by: and, or, xor, eq, ge, gt, le, lt, ne
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  z = cadaEmptyEval(x,y);
  return
end
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;

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
end

% ----------------------------Function Sizing---------------------------- %
if (xMrow == yMrow && xNcol == yNcol)
  FMrow = yMrow; FNcol = yNcol;
elseif (xMrow == 1 && xNcol == 1)
  FMrow = yMrow; FNcol = yNcol;
elseif (yMrow == 1 && yNcol == 1)
  FMrow = xMrow; FNcol = xNcol;
elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
  % In printing run of a FOR loop - Sizes could be off due to IF statement
  if xMrow <= yMrow && xNcol <= yNcol
    y = cadaRemoveRowsCols(y,[xMrow xNcol]);
    yMrow = xMrow; yNcol = xNcol;
  elseif yMrow <= xMrow && yNcol <= xNcol
    x = cadaRemoveRowsCols(x,[yMrow yNcol]);
    xMrow = yMrow; xNcol = yNcol;
  else
    error('Inputs are not of compatible sizes');
  end
  FMrow = xMrow; FNcol = xNcol;
else
  error('Inputs are not of compatible sizes');
end

% ----------------------Build Function Properties--------------------------
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[],'logical',[]);
if isinf(FMrow)
  FMrow = 1; xMrow = 1; yMrow = 1;
elseif isinf(FNcol)
  FNcol = 1; xNcol = 1; yNcol = 1;
end

%----------------------Function Numeric Values (if exist)-----------------%
if ~isempty(x.func.value) && ~isempty(y.func.value)
  callerfun = str2func(callerstr);
  z.func.value = callerfun(x.func.value,y.func.value);
else
  if ~isempty(x.func.value)
    xtemp = logical(x.func.value);
  elseif ~isempty(x.func.zerolocs)
    xtemp = true(xMrow,xNcol);
    xtemp(x.func.zerolocs) = false;
  else
    xtemp = true(xMrow,xNcol);
  end
  if ~isempty(y.func.value)
    ytemp = logical(y.func.value);
  elseif ~isempty(x.func.zerolocs)
    ytemp = true(yMrow,yNcol);
    ytemp(x.func.zerolocs) = false;
  else
    ytemp = true(yMrow,yNcol);
  end
  switch callerstr
    case 'and'
      % Element i of z is known to be zero if either element i of x or y
      % is known to be zero.
      ztemp = and(xtemp,ytemp);
      z.func.zerolocs = find(~ztemp(:));
    case 'or'
      % Element i of z is known to be zero if element i of x AND element i
      % of y is known to be zero
      ztemp = or(xtemp,ytemp);
      z.func.zerolocs = find(~ztemp(:));
    case 'xor'
      % Only way we can know for certain that an element of z is zero is
      % that we need the actual numeric values of either x or y together
      % with a zero loc of the other.
      if ~isempty(x.func.value) || ~isempty(y.func.value)
        ztemp = xor(xtemp,ytemp);
      else
        ztemp = true(FMrow,FNcol);
      end
    case 'eq'
      % Only way we can that two elements are not equal is that if we know
      % the numeric values of one and the zero values of the other.
      ztemp = true(FMrow,FNcol);
      if ~isempty(x.func.value)
        ztemp(~ytemp & xtemp) = false;
      elseif ~isempty(y.func.value)
        ztemp(~xtemp & ytemp) = false;
      end
    case 'ne'
      % If element i of x and y are both zero then we know that element i
      % of z is also zero
      ztemp = true(FMrow,FNcol);
      ztemp(~xtemp & ~ytemp) = false;
    case {'gt','ge'}
      ztemp = true(FMrow,FNcol);
      if ~isempty(x.func.value)
        ztemp(~ytemp & x.func.value < 0) = false;
      elseif ~isempty(y.func.value)
        ztemp(~xtemp & y.func.value > 0) = false;
      end
    case {'lt','le'}
      ztemp = true(FMrow,FNcol);
      if ~isempty(x.func.value)
        ztemp(x.func.value > 0 & ~ytemp) = false;
      elseif ~isempty(y.func.value)
        ztemp(y.func.value < 0 & ~xtemp) = false;
      end
  end
  z.func.zerolocs = find(~ztemp(:));
  if length(z.func.zerolocs) == FMrow*FNcol
    z.func.zerolocs = []; z.func.value = false(FMrow,FNcol);
  end
end

% ----------------------------Function Printing ------------------------- %
if PFLAG == 1
  fprintf(fid,[indent,funcstr,' = ',callerstr,'(',x.func.name,',',y.func.name,');\n']);
end


% ------------------------Build Derivative Properties----------------------
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
