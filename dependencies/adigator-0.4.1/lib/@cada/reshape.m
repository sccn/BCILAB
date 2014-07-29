function y = reshape(x,varargin)
% CADA overloaded version of function RESHAPE.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
if ADIGATOR.FORINFO.FLAG
  IncreaseForReshapeCount();
  if ADIGATOR.RUNFLAG == 2
    [y,flag,x] = ForReshape(x,varargin);
    if flag; return; end
  end
end
if ADIGATOR.EMPTYFLAG
  % Need to send to cadaEmptyEval, build the string to call it
  InputStr    = cell(1,nargin);
  InputStr{1} = 'x,';
  for Icount = 2:nargin
    InputStr{Icount} = sprintf('varargin{%1.0f},',Icount-1);
  end
  InputStr = cell2mat(InputStr);
  ForOpCall = ['cadaEmptyEval(',InputStr(1:end-1),');'];
  y = eval(ForOpCall);
  return
end
% ----------------------------Determine Inputs----------------------------%
varargSize = length(varargin);
if varargSize == 0 ;
  error('Requires at least 2 inputs.');
elseif varargSize == 1
  if isa(varargin{1},'cada')
    if ~isempty(varargin{1}.func.value)
      ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
      if length(varargin{1}.func.value) > 2
        error('cannot RESHAPE more than 2 dimensions')
      elseif length(varargin{1}.func.value) == 1
        error('Size vector must have at least two elements.');
      end
      FMrow = varargin{1}.func.value(1);
      FNcol = varargin{1}.func.value(2);
      RepStr = varargin{1}.func.name;
    else
      error('cannot RESHAPE a strictly symbolic dimension')
    end
  elseif length(varargin{1}) == 2
    FMrow  = varargin{1}(1);
    FNcol  = varargin{1}(2);
    RepStr = sprintf('%1.0f,%1.0f',FMrow,FNcol);
  else
    error('cannot RESHAPE more than 2 dimensions')
  end
elseif varargSize == 2
  FMrow = varargin{1};
  if isa(FMrow,'cada')
    if ~isempty(FMrow.func.value)
      ADIGATOR.VARINFO.LASTOCC(FMrow.id,1) = ADIGATOR.VARINFO.COUNT;
      RowStr = FMrow.func.name;
      FMrow  = FMrow.func.value;
    else
      error('cannot RESHAPE a strictly symbolic dimension')
    end
  else
    RowStr = sprintf('%1.0f',FMrow);
  end
  FNcol = varargin{2};
  if isa(FNcol,'cada')
    if ~isempty(FNcol.func.value)
      ADIGATOR.VARINFO.LASTOCC(FNcol.id,1) = ADIGATOR.VARINFO.COUNT;
      ColStr = FNcol.func.name;
      FNcol  = FNcol.func.value;
    else
      error('cannot RESHAPE a strictly symbolic dimension')
    end
  else
    ColStr = sprintf('%1.0f',FNcol);
  end
  RepStr = ['',RowStr,',',ColStr,''];
else
  error('Too many input arguments.');
end
if FMrow*FNcol ~= x.func.size(1)*x.func.size(2)
  error('To RESHAPE the number of elements must not change.')
end

% --------------------------Parse X-------------------------------------- %
if isnumeric(x)
    [xMrow,xNcol] = size(x);
    xtemp.id = [];
    xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
      'value',x);
    if PFLAG
      xtemp.func.name = cadamatprint(x);
    end
    xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    x = xtemp;
elseif ~isa(x,'cada')
    y = repmat(x,FMrow,FNcol);
    return
else
  xMrow = x.func.size(1);
  xNcol = x.func.size(2);
end

% --------------------Build Function Properties---------------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow FNcol],...
  'zerolocs',x.func.zerolocs,'value',[]);
% Vectorized stuff
if isinf(FMrow)  && ...
    ((isinf(xMrow) && xNcol == 1) || (xMrow == 1 && isinf(xNcol)))
  FMrow = 1;
elseif isinf(FNcol)  && ...
    ((isinf(xMrow) && xNcol == 1) || (xMrow == 1 && isinf(xNcol)))
  FNcol = 1;
elseif isinf(FMrow) || isinf(FNcol)
  error(['Can only reshape vectorized input if input is of size N x 1 ',...
    'or 1 x N, where N is the vectorized dimension']);
end
if ~isempty(x.func.value)
  y.func.value = reshape(x.func.value,FMrow,FNcol);
end
if isfield(x.func,'logical')
  y.func.logical = [];
end

% ------------------Build Derivative Properties-----------------------%
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name   = derivstr;
    y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
    % ----------------Derivative Printing-------------------------%
    if DPFLAG
      fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
    end
  end
end

% -----------------------Function Printing ---------------------------%
if PFLAG
  fprintf(fid,[indent,funcstr,' = reshape(',x.func.name,',',RepStr,');\n']);
end
ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG
  AssignForReshapeData(y,x);
end
return
end

function IncreaseForReshapeCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.RESHAPE = ...
  ADIGATORFORDATA(INNERLOC).COUNT.RESHAPE + 1;
return
end

function AssignForReshapeData(y,x)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Rcount    = ADIGATORFORDATA(INNERLOC).COUNT.RESHAPE;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;

x.func.size(isinf(x.func.size)) = 1;
y.func.size(isinf(y.func.size)) = 1;
% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES =...
    [y.func.size.';x.func.size.'];
else
  ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES(:,ITERCOUNT) = ...
    [y.func.size.';x.func.size.'];
end

% Variable OverMapping
if ~isa(x,'cada'); x = class(x,'cada'); end
if isempty(ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS)
  % First Call
  ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{1} = y;
  if ~isempty(x.id) && ADIGATOR.VARINFO.NAMELOCS(x.id)
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = x;
    end
  else
    ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = x;
  end
else
  yOver = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{1};
  ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{1} = cadaUnionVars(y,yOver);
  xOver = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2};
  if isa(xOver,'cada')
    ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = cadaUnionVars(x,xOver);
  else
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    if xOverLoc1 && xOver ~= xOverLoc1
      ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2} = xOverLoc1;
    end
  end
end
return
end

function [y,flag,x] = ForReshape(x,Sinputs)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Rcount    = ADIGATORFORDATA(INNERLOC).COUNT.RESHAPE;

xOver = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{2};
% Check that X is overmapped properly
x = cadaPrintReMap(x,xOver,x.id);

if isempty(ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES)
  flag = 0; y = []; return
else
  flag = 1;
end
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NUMvod    = ADIGATOR.NVAROFDIFF;
NDstr     = sprintf('%1.0d',ADIGATOR.DERNUMBER);

%--------------------Parse Inputs---------------------------------------- %
if length(Sinputs) == 1
  if isa(Sinputs{1},'cada')
    if Sinputs{1}.func.size(1)*Sinputs{1}.func.size(2) == 1
      resMstr = Sinputs{1}.func.name;
      ResStr = [resMstr,',',resMstr];
    else
      ResStr = Sinputs{1}.func.name;
    end
    ADIGATOR.VARINFO.LASTOCC(Sinputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{1})
    if length(Sinputs{1}) == 1
      resMrow = Sinputs{1};
      resNcol = resMrow;
    elseif length(Sinputs{1}) == 2
      resMrow = Sinputs{1}(1);
      resNcol = Sinputs{1}(2);
    end
    ResStr = sprintf('%1.0d,%1.0d',resMrow,resNcol);
  end
elseif length(Sinputs) == 2
  if isa(Sinputs{1},'cada')
    resMstr = Sinputs{1}.func.name;
    ADIGATOR.VARINFO.LASTOCC(Sinputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{1})
    resMstr = sprintf('%1.0d',Sinputs{1});
  end
  if isa(Sinputs{2},'cada')
    resNstr = Sinputs{2}.func.name;
    ADIGATOR.VARINFO.LASTOCC(Sinputs{2}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{2})
    resNstr = sprintf('%1.0d',Sinputs{2});
  end
  ResStr =[resMstr,',',resNstr];
end

% Get OverMapped Variables
yOver = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).VARS{1};


% Build Y
y    = yOver;
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func.name      = funcstr;


% ----------------Print Out Derivatives---------------------------------- %
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
      IndName = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).INDICES{Vcount,1};
      DepFlag = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).INDICES{Vcount,3}(1);
      SpFlag  = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).INDICES{Vcount,3}(2);
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

% ----------------------Print Out Function------------------------------- %
TF1 = ['cada',NDstr,'tempf1'];
% Get Our Reference off of X first
RowInds = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{1,1};
ColInds = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{2,1};
if ~isempty(RowInds) && ~isempty(ColInds)
  RowDepFlag = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{1,3}(1);
  ColDepFlag = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{2,3}(1);
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
  xStr = [x.func.name,'(1:',RowRef,',1:',ColRef,')'];
elseif ~isempty(RowInds)
  RowDepFlag = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{1,3}(1);
  if RowDepFlag
    RowRef = [RowInds,'(',CountName,')'];
  else
    RowRef = RowInds;
  end
  if x.func.size(2) == 1
    xStr = [x.func.name,'(1:',RowRef,')'];
  else
    xStr = [x.func.name,'(1:',RowRef,',:)'];
  end
elseif ~isempty(ColInds)
  ColDepFlag = ADIGATORFORDATA(INNERLOC).RESHAPE(Rcount).SIZES{2,3}(1);
  if ColDepFlag
    ColRef = [ColInds,'(',CountName,')'];
  else
    ColRef = ColInds;
  end
  if x.func.size(1) == 1
    xStr = [x.func.name,'(1:',ColRef,')'];
  else
    xStr = [x.func.name,'(:,1:',ColRef,')'];
  end
else
  xStr = x.func.name;
end
% Print out the Temp RESHAPE
fprintf(fid,[indent,TF1,' = reshape(',xStr,',',ResStr,');\n']);
% Print out Y
fprintf(fid,[indent,funcstr,' = zeros(%1.0d,%1.0d);\n'],y.func.size(1),y.func.size(2));
fprintf(fid,[indent,funcstr,'(1:size(',TF1,',1),1:size(',TF1,',2)) = ',TF1,';\n']);

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
  y = class(y,'cada');
end
return
end
