function y = sparse(varargin)
% CADA overloaded version of MATLAB function SPARSE.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
indent  = ADIGATOR.PRINT.INDENT;
PFLAG   = ADIGATOR.PRINT.FLAG;
if ADIGATOR.FORINFO.FLAG && nargin > 1
  IncreaseForSparseCount();
  if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    y = ForSparse(varargin);
    return
  end
end
if ADIGATOR.EMPTYFLAG
  if nargin > 1
    % Need to construct a call to the for op holder
    InputStr = cell(1,nargin);
    for Icount = 1:nargin
      InputStr{Icount} = sprintf('varargin{%1.0d},',Icount);
    end
    InputStr = cell2mat(InputStr);
    ForOpStr = ['cadaEmptyEval(',InputStr(1:end-1),');'];
    y = eval(ForOpStr);
  else
    y = cadaEmptyEval(varargin{1});
  end
  return
end
% ------------------------Parse Inputs----------------------------------- %
if nargin == 1
  % y = sparse(x);
  x = varargin{1};
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,DPFLAG] = cadafuncname();
  y.func = x.func;
  y.func.name = funcstr;
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
      if DPFLAG
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
    end
  end
  if PFLAG
    fprintf(fid,[indent,funcstr,' = sparse(',x.func.name,');\n']);
  end
  ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');
  return
elseif nargin == 3 || nargin == 5 || nargin == 6
  yrows = varargin{1};  ycols = varargin{2};  x = varargin{3};
  if isa(yrows,'cada')
    if ~isempty(yrows.func.value)
      ADIGATOR.VARINFO.LASTOCC(yrows.id,1) = ADIGATOR.VARINFO.COUNT;
      RowStr = yrows.func.name;
      yrows  = yrows.func.value;
    else
      error('cannot index with a strictly symbolic variable')
    end
  else
    if PFLAG; RowStr = cadaindprint(yrows); else RowStr = []; end
  end
  if isa(ycols,'cada')
    if ~isempty(ycols.func.value)
      ADIGATOR.VARINFO.LASTOCC(ycols.id,1) = ADIGATOR.VARINFO.COUNT;
      ColStr = ycols.func.name;   ycols = ycols.func.value;
    else
      error('cannot index with a strictly symbolic variable')
    end
  else
    if PFLAG; ColStr = cadaindprint(ycols); else ColStr = []; end
  end
  if isa(x,'cada')
    xStr = x.func.name;
  else
    if PFLAG; xStr = cadamatprint(x); else xStr = []; end
    xtemp.id = [];
    xtemp.func = struct('name',xStr,'size',size(x),'zerolocs',[],...
      'value',x);
    xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    x = xtemp;
  end
  if nargin == 3
    FMrow = max(yrows);   FNcol = max(ycols);
    OutStr = ['sparse(',RowStr,',',ColStr,',',xStr,')'];
  else
    FMrow = varargin{4};  FNcol = varargin{5};
    if isa(FMrow,'cada')
      if ~isempty(FMrow.func.value)
        ADIGATOR.VARINFO.LASTOCC(FMrow.id,1) = ADIGATOR.VARINFO.COUNT;
        FMrowStr  = FMrow.func.name;
        FMrow     = FMrow.func.value;
      else
        error('cannot index with a strictly symbolic variable')
      end
    else
      FMrowStr = sprintf('%1.0f',FMrow);
    end
    if isa(FNcol,'cada')
      if ~isempty(FNcol.func.value)
        ADIGATOR.VARINFO.LASTOCC(FNcol.id,1) = ADIGATOR.VARINFO.COUNT;
        FNcolStr = FNcol.func.name;
        FNcol = FNcol.func.value;
      else
        error('cannot index with a strictly symbolic variable')
      end
    else
      FNcolStr = sprintf('%1.0f',FNcol);
    end
    if nargin == 5
      OutStr = ['sparse(',RowStr,',',ColStr,',',xStr,',',FMrowStr,',',FNcolStr,')'];
    else
      NZmax = varargin{6};
      if isa(NZmax,'cada')
        if ~isempty(NZmax.func.value)
          ADIGATOR.VARINFO.LASTOCC(NZmax.id,1) = ADIGATOR.VARINFO.COUNT;
          NZmaxStr = NZmax.func.name;
        else
          error('cannot index with a strictly symbolic variable')
        end
      else
        NZmaxStr = sprintf('%1.0f',NZmax);
      end
      OutStr = ['sparse(',RowStr,',',ColStr,',',xStr,',',FMrowStr,',',FNcolStr,',',NZmaxStr,')'];
    end
  end
else
  error('invalid number of inputs')
end

% ------------------------Get Y------------------------------------------ %
yrows = yrows(:); ycols = ycols(:);
nsubs = length(yrows);
if nsubs ~= length(ycols) || nsubs ~= x.func.size(1)*x.func.size(2)
  error('Vectors must be the same lengths.')
end
isubs = sub2ind([FMrow,FNcol],yrows,ycols);
if nsubs ~= length(unique(isubs))
  error('Currently Cannot do SPARSE command when indices are not unique')
end
if size(isubs,2) > 1
  isubs = isubs';
end

% -----------------Build Function Properties-----------------------------
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);

% --------------------Function Numerics/Sparsity---------------------------
if ~isempty(x.func.value)
  y.func.value = sparse(yrows,ycols,x.func.value,FMrow,FNcol);
else
  ytemp           = false(FMrow,FNcol);
  xtemp = true(length(isubs),1);
  if ~isempty(x.func.zerolocs)
    xtemp(x.func.zerolocs) = false;
  end
  ytemp(isubs)    = xtemp;
  y.func.zerolocs = find(~ytemp(:));
  if length(y.func.zerolocs) == FMrow*FNcol
    y.func.zerolocs = [];
  end
end



% --------------------Build Derivative Properties------------------------
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    xrows = x.deriv(Vcount).nzlocs(:,1);
    xcols = x.deriv(Vcount).nzlocs(:,2);
    xrows = isubs(xrows);
    dy = sparse(xrows,xcols,1:nnz(xrows),FMrow*FNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
    [yrows,ycols,yind] = find(dy);
    if size(yrows,2) > 1
      yrows = yrows';
      ycols = ycols';
    end
    y.deriv(Vcount).nzlocs = [yrows,ycols];
    if DPFLAG
      Dind1 = cadaindprint(yind);
      fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
    end
  end
end

if PFLAG
  fprintf(fid,[funcstr,' = ',OutStr,';\n']);
end

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG && nargin > 1
  AssignForSparseData(y,x,isubs);
end
return
end

function IncreaseForSparseCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.SPARSE =...
  ADIGATORFORDATA(INNERLOC).COUNT.SPARSE+1;
return
end

function AssignForSparseData(y,x,inds)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
SPcount   = ADIGATORFORDATA(INNERLOC).COUNT.SPARSE;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;

% Indice Assignment
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).INDICES = inds(:);
else
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).INDICES(1:length(inds),ITERCOUNT) = inds;
end

% Variable OverMapping
if ~isa(x,'cada')
  x = class(x,'cada');
end

if isempty(ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS)
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{1} = y;
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).SIZES = [y.func.size.';x.func.size.'];
  if ~isempty(x.id) && ADIGATOR.VARINFO.NAMELOCS(x.id)
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = x;
    end
  else
    ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = x;
  end
else
  yOver = ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{1};
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{1} = cadaUnionVars(y,yOver);
  ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).SIZES(:,end) = [y.func.size.';x.func.size.'];
  xOver = ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2};
  if isa(xOver,'cada')
    ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = cadaUnionVars(x,xOver);
  else
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    if xOverLoc1 && xOver ~= xOverLoc1
      ADIGATORFORDATA(INNERLOC).SPARSE(SPcount).VARS{2} = xOverLoc1;
    end
  end
end
return
end

function y = ForSparse(Inputs)
global ADIGATOR ADIGATORFORDATA
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0d',ADIGATOR.DERNUMBER);
NUMvod    = ADIGATOR.NVAROFDIFF;
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Scount    = ADIGATORFORDATA(INNERLOC).COUNT.SPARSE;
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
% ----------------------Parse Inputs------------------------------------- %
NUMinputs = length(Inputs);
if NUMinputs == 3 || NUMinputs == 5 || NUMinputs == 6
  yrows = Inputs{1};
  ycols = Inputs{2};
  x     = Inputs{3};
  if isa(yrows,'cada')
    ADIGATOR.VARINFO.LASTOCC(yrows.id,1) = ADIGATOR.VARINFO.COUNT;
    RowStr = yrows.func.name;
  else
    RowStr = cadaindprint(yrows);
  end
  if isa(ycols,'cada')
    ADIGATOR.VARINFO.LASTOCC(ycols.id,1) = ADIGATOR.VARINFO.COUNT;
    ColStr = ycols.func.name;
  else
    ColStr = cadaindprint(ycols);
  end
  if isa(x,'cada')
    xStr = x.func.name;
  else
    xtemp.id = [];
    xStr = cadamatprint(x);
    xtemp.func = struct('name',xStr,'size',size(x),'zerolocs',[],...
      'value',x);
    xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    x = xtemp;
  end
  if NUMinputs > 3
    FMrow = Inputs{4};
    FNcol = Inputs{5};
    if isa(FMrow,'cada')
      ADIGATOR.VARINFO.LASTOCC(FMrow.id,1) = ADIGATOR.VARINFO.COUNT;
    end
    if isa(FNcol,'cada')
      ADIGATOR.VARINFO.LASTOCC(FNcol.id,1) = ADIGATOR.VARINFO.COUNT;
    end
    if NUMinputs == 6
      NZmax = Inputs{6};
      if isa(NZmax,'cada')
        ADIGATOR.VARINFO.LASTOCC(NZmax.id,1) = ADIGATOR.VARINFO.COUNT;
      else
      end
    end
  end
else
  error('invalid number of inputs')
end

% Get OverMapped X and Y
yOver = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).VARS{1};
xOver = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).VARS{2};

% Check and Make sure X is properly OverMapped.
x     = cadaPrintReMap(x,xOver,x.id);

% Built y function
y    = yOver;
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func.name      = funcstr;



% -------------------Print out Derivatives------------------------------- %
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
      IndName = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).INDICES{Vcount,1};
      DepFlag = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).INDICES{Vcount,3}(1);
      SpFlag  = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).INDICES{Vcount,3}(2);
      if DepFlag
        % Indices are dependent upon this loop
        IndRef = [IndName,'(:,',CountName,')'];
      else
        % Indices are independent of this loop
        IndRef = IndName;
      end
      if SpFlag
        % Indices are sparse
        fprintf(fid,[indent,TD1,'(logical(',IndRef,')) = ',...
          x.deriv(Vcount).name,'(nonzeros(',IndRef,'));\n']);
      else
        % Indices are non-sparse
        fprintf(fid,[indent,TD1,' = ',...
          x.deriv(Vcount).name,'(',IndRef,');\n']);
      end
      fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
    end
  end
end

%---------------------Print Out Function--------------------------------- %
RowInds = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{1,1};
ColInds = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{2,1};
% If X is changing Row or Cols, need to reference off of X and the
% Rows/Cols inputs
if ~isempty(RowInds) && ~isempty(ColInds)
  % Input X is changing both Row and Column Sizes
  RowDepFlag = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{1,3}(1);
  ColDepFlag = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{2,3}(1);
  if RowDepFlag
    RowRef = [RowInds,'(',CountName,')'];
  else
    RowRef = RowInds;
  end
  if ColDepFlag
    ColRef = [ColInds,'(',CountName,')'];
  else
    ColRef = ColInds;
  end
  xStr = [xStr,'(1:',RowRef,',1:',ColRef,')'];
  RowStr = [RowStr,'(1:',RowRef,'*',ColRef,')'];
  ColStr = [ColStr,'(1:',RowRef,'*',ColRef,')'];
elseif ~isempty(RowInds)
  % Input X is changing Row Sizes
  RowDepFlag = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{1,3}(1);
  if RowDepFlag
    RowRef = [RowInds,'(',CountName,')'];
  else
    RowRef = RowInds;
  end
  if x.func.size(2) == 1
    xStr   = [xStr,'(1:',RowRef,')'];
    RowStr = [RowStr,'(1:',RowRef,')'];
    ColStr = [ColStr,'(1:',RowRef,')'];
  else
    xStr   = [xStr,'(1:',RowRef,',:)'];
    RowStr = sprintf([RowStr,'(1:',RowRef,'*%1.0d)'],x.func.size(2));
    ColStr = sprintf([ColStr,'(1:',RowRef,'*%1.0d)'],x.func.size(2));
  end
elseif ~isempty(ColInds)
  % Input X is changing Column Sizes
  ColDepFlag = ADIGATORFORDATA(INNERLOC).SPARSE(Scount).SIZES{2,3}(1);
  if ColDepFlag
    ColRef = [ColInds,'(',CountName,')'];
  else
    ColRef = ColInds;
  end
  if x.func.size(1) == 1
    xStr   = [xStr,'(1:',ColRef,')'];
    RowStr = [RowStr,'(1:',ColRef,')'];
    ColStr = [ColStr,'(1:',ColRef,')'];
  else
    xStr   = [xStr,'(:,1:',ColRef,')'];
    RowStr = sprintf([RowStr,'(1:',ColRef,'*%1.0d)'],x.func.size(1));
    ColStr = sprintf([ColStr,'(1:',ColRef,'*%1.0d)'],x.func.size(1));
  end
end
% Print out the function
fprintf(fid,[indent,funcstr,' = sparse(',RowStr,',',ColStr,',',xStr,...
  ',%1.0d,%1.0d);\n'],y.func.size(1),y.func.size(2));

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
y = class(y,'cada');
end
return
end