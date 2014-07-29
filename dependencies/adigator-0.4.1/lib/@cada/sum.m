function y = sum(x,varargin)
% CADA overloaded version of function SUM
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  if nargin == 1
    y = cadaEmptyEval(x);
  else
    y = cadaEmptyEval(x,varargin{1});
  end
  return
end
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% ----------------------------Parse Inputs------------------------------- %
if nargin == 1
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  if xMrow == 1; Dim = 2; else Dim = 1; end
elseif nargin == 2
  if isa(x,'cada')
    xMrow = x.func.size(1); xNcol = x.func.size(2);
  elseif isnumeric(x)
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
  else
    error('??? Cannot sum a non-numeric object.')
  end
  Dim = varargin{1};
  if isa(Dim,'cada')
    if ~isempty(Dim.func.value)
      ADIGATOR.VARINFO.LASTOCC(Dim.id,1) = ...
        ADIGATOR.VARINFO.COUNT;
      Dim = Dim.func.value;
    else
      error('??? cannot sum across a purely symbolic dimension')
    end
  elseif ~isnumeric(varargin{1})
    error('??? Invalid input argument.')
  end
  
else
  error('Too many input arguments')
end
if Dim == 1
  yMrow = 1; yNcol = xNcol;
else
  yMrow = xMrow; yNcol = 1;
end
% -----------------------Build Y Function-------------------------------- %
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[yMrow yNcol],'zerolocs',[],'value',[]);
if isinf(xMrow) && Dim == 2
  xMrow = 1; xvec = 1;
elseif isinf(xNcol) && Dim == 1
  xNcol = 1; xvec = 2;
elseif isinf(xMrow) || isinf(xNcol)
  error('Cannot sum over vectorized dimension')
else
  xvec = 0;
end


% Function Numerics/Sparsity
if ~isempty(x.func.value)
  % Y is numeric
  y.func.value = sum(x.func.value,Dim);
elseif ~isempty(x.func.zerolocs)
  % X is sparse
  xtemp = true(xMrow,xNcol);
  xtemp(x.func.zerolocs) = false;
  ytemp = sum(xtemp,Dim);
  y.func.zerolocs = find(~ytemp(:));
end

% --------------------------Build Derivative----------------------------- %
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if cadaCheckForDerivs(x)
  % x has derivatives - going to create a dummy variable and call
  % cadamtimesderiv
  TF2 = ['cada',NDstr,'tf2'];
  if Dim == 1
    % A is 1 by xMrow
    A.func.size = [1 xMrow];
    % - sum(x,1) is same as A*X
  else
    % A is xNcol by 1
    A.func.size = [xNcol 1];
    % - sum(x,2) is same as X*A
  end
  A.func.name = TF2;
  if DPFLAG
    fprintf(fid,[indent,TF2,' = ones(%1.0f,%1.0f);\n'],...
      A.func.size(1),A.func.size(2));
  end
  Atemp = ones(A.func.size);
  A.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  xtemp = ones(xMrow,xNcol);
end
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    if Dim == 1
      % y = A*X
      if xvec
        nzlocs = cadamtimesderivvec(A,x,Atemp,xtemp,Vcount,derivstr,DPFLAG,xvec);
      else
        nzlocs = cadamtimesderiv(A,x,Atemp,xtemp,Vcount,derivstr,DPFLAG,'mtimes');
      end
    else
      % y = X*A
      if xvec
        nzlocs = cadamtimesderivvec(x,A,xtemp,Atemp,Vcount,derivstr,DPFLAG,xvec);
      else
        nzlocs = cadamtimesderiv(x,A,xtemp,Atemp,Vcount,derivstr,DPFLAG,'mtimes');
      end
    end
    y.deriv(Vcount).nzlocs = nzlocs;
    y.deriv(Vcount).name = derivstr;
  end
end


% ---------------------Print Out Function-------------------------------- %
if PFLAG
  if nargin == 1
    fprintf(fid,[indent,funcstr,' = sum(',x.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = sum(',x.func.name,',%1.0d);\n'],Dim);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;