function y = diag(x,varargin)
% CADA overloaded version of function DIAG
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;

if ADIGATOR.FORINFO.FLAG
  IncreaseForOtherCount();
end
if ADIGATOR.EMPTYFLAG
  if nargin == 1
    y = cadaEmptyEval(x);
  else
    y = cadaEmptyEval(x,varargin{1});
  end
  return
end
% ---------------------------- Parse Inputs ------------------------------%
if nargin == 2
  if isa(x,'cada')
    xMrow = x.func.size(1);
    xNcol = x.func.size(2);
  else
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
  if isa(varargin{1},'cada')
    K = varargin{1}.func.value;
    if isempty(K)
      error('cannot assign to strictly symbolic diagonal')
    end
    Kstr = varargin{1}.func.name;
  else
    K    = varargin{1};
    Kstr = sprintf('%1.0f',K);
  end
else
  xMrow = x.func.size(1);
  xNcol = x.func.size(2);
  K    = 0;
  Kstr = '0';
end

if xMrow == 1 || xNcol == 1
  % --------------------------------------------------------------------- %
  %                          Assignment Case                              %
  % --------------------------------------------------------------------- %
  if xMrow == 1; N = xNcol; else N = xMrow; end
  if isinf(N)
    error('Cannot use diag on vectorized object')
  end
  
  % -------------------- Build Function Properties ---------------------- %
  FMrow = abs(K) + N; FNcol = FMrow;
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,DPFLAG] = cadafuncname();
  y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
    'value',[]);
  
  if ~isempty(x.func.value)
    xtemp = x.func.value;
    y.func.value = diag(xtemp,K);
  else
    xtemp = ones(N,1);
    if ~isempty(x.func.zerolocs)
      xtemp(x.func.zerolocs) = 0;
    end
    ytemp = diag(xtemp,K);
    y.func.zerolocs = find(~logical(ytemp(:)));
  end
  
  yref    = zeros(FMrow,FNcol);
  yref(:) = 1:FMrow*FNcol;
  xref    = diag(yref,K);
  % -------------------- Build Derivative Properties -------------------- %
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  for Vcount = 1:NUMvod
    dxind = x.deriv(Vcount).nzlocs;
    if ~isempty(dxind)
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name   = derivstr;
      yrows = xref(dxind(:,1)); ycols = dxind(:,2);
      y.deriv(Vcount).nzlocs = [yrows ycols];
      if DPFLAG
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
    end
  end
else
  % --------------------------------------------------------------------- %
  %                           Reference Case                              %
  % --------------------------------------------------------------------- %
  xtemp = ones(xMrow,xNcol);
  ytemp = diag(xtemp,K);
  [FMrow FNcol] = size(ytemp);
  
  % -------------------- Build Function Properties ---------------------- %
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,DPFLAG] = cadafuncname();
  y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
    'value',[]);
  
  if ~isempty(x.func.value)
    xtemp = x.func.value;
    y.func.value = diag(xtemp,K);
  else
    xtemp = ones(xMrow,xNcol);
    if ~isempty(x.func.zerolocs)
      xtemp(x.func.zerolocs) = 0;
    end
    ytemp = diag(xtemp,K);
    y.func.zerolocs = find(~logical(ytemp(:)));
  end
  
  xref    = zeros(xMrow,xNcol);
  xref(:) = 1:xMrow*xNcol;
  yref    = diag(xref,K);
  % -------------------- Build Derivative Properties -------------------- %
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  for Vcount = 1:NUMvod
    dxind = x.deriv(Vcount).nzlocs;
    if ~isempty(dxind)
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name   = derivstr;
      dx = sparse(dxind(:,1),dxind(:,2),1:size(dxind,1),xMrow*xNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
      dy = dx(yref,:);
      [yrows ycols yind] = find(dy);
      if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
      y.deriv(Vcount).nzlocs = [yrows ycols];
      if DPFLAG
        Dind1 = cadaindprint(yind(:));
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
      end
    end
  end
end

% -----------------------Function Printing ---------------------------%
if PFLAG
  fprintf(fid,[indent,funcstr,' = diag(',x.func.name,',',Kstr,');\n']);
end
ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG
  ForOtherData([xMrow xNcol FMrow FNcol K])
end
return
end

function IncreaseForOtherCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.OTHER = ADIGATORFORDATA(INNERLOC).COUNT.OTHER+1;
return
end

function ForOtherData(data)
global ADIGATOR ADIGATORFORDATA
% Store:
%       1. What case the logical statement lies in
%       2. What the size of the logical statement is
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.OTHER;
oldData = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA;
if isempty(oldData)
  ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA = data;
elseif any(oldData ~= data)
  error('Currenlty diag is not coded such that the inputs may change size or diagonal within a loop')
end
return
end
