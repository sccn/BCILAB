function y = cadacreatearray(callerstr,Inputs)
% Called by zeros, ones, eye, true, false
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
  if ADIGATOR.RUNFLAG == 2;
    y = ForZeros(Inputs,callerstr);
    return
  end
end

if length(Inputs) == 1
  if ADIGATOR.EMPTYFLAG
    y = cadaEmptyEval(Inputs{1});
    return
  end
  if ~isempty(Inputs{1}.func.value)
    ADIGATOR.VARINFO.LASTOCC(Inputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
    if Inputs{1}.func.size(1)*Inputs{1}.func.size(2) == 1
      FMrow = Inputs{1}.func.value(1);
      FNcol = FMrow;
    elseif Inputs{1}.func.size(1)*Inputs{1}.func.size(2) == 2
      FMrow = Inputs{1}.func.value(1);
      FNcol = Inputs{1}.func.value(2);
    else
      error('can only be used with 2 dimensions')
    end
    SizeStr = Inputs{1}.func.name;
  else
    error('cannot pre-allocate a stricly Symbolic size')
  end
elseif length(Inputs) == 2
  if ADIGATOR.EMPTYFLAG
    y = cadaEmptyEval(Inputs{1},Inputs{2});
    return
  end
  if isa(Inputs{1},'cada')
    if ~isempty(Inputs{1}.func.value)
      ADIGATOR.VARINFO.LASTOCC(Inputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
      FMrow = Inputs{1}.func.value(1);
      FMstr = Inputs{1}.func.name;
    else
      error('cannot pre-allocate a stricly Symbolic size')
    end
  elseif isnumeric(Inputs{1})
    FMrow = Inputs{1};
    FMstr = sprintf('%1.0f',FMrow);
  else
    error('invalid input')
  end
  if isa(Inputs{2},'cada')
    if ~isempty(Inputs{2}.func.value)
      ADIGATOR.VARINFO.LASTOCC(Inputs{2}.id,1) = ADIGATOR.VARINFO.COUNT;
      FNcol = Inputs{2}.func.value(1);
      FNstr = Inputs{2}.func.name;
    else
      error('cannot pre-allocate a stricly Symbolic size')
    end
  elseif isnumeric(Inputs{2})
    FNcol = Inputs{2};
    FNstr = sprintf('%1.0f',FNcol);
  else
    error('invalid input')
  end
  SizeStr = [FMstr,',',FNstr];
else
  error('can only be used with 2 dimensions')
end

% --------------------Build Function Properties-----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);
if isinf(FMrow) && isinf(FNcol)
  error('At least one dimension must be non-vectorized')
elseif isinf(FMrow)
  FMrow = 1; yvec = 1;
elseif isinf(FNcol)
  FNcol = 1; yvec = 2;
else
  yvec = 0;
end

y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
switch callerstr
  case 'zeros'
    y.func.value = zeros(FMrow,FNcol);
  case 'ones'
    y.func.value = ones(FMrow,FNcol);
  case 'eye'
    if ~yvec
      y.func.value = eye(FMrow,FNcol);
    end
  case 'true'
    y.func.value = true(FMrow,FNcol);
  case 'false'
    y.func.value = false(FMrow,FNcol);
end

if PFLAG == 1
  fprintf(fid,[indent,funcstr,' = ',callerstr,'(',SizeStr,');\n']);
end

ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;

if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 1
  ForOtherData(y.func.size);
end
y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
return
end

function IncreaseForOtherCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.OTHER = ADIGATORFORDATA(INNERLOC).COUNT.OTHER+1;
return
end

function ForOtherData(y)
% Storing Maximum Size.
global ADIGATOR ADIGATORFORDATA
INNERLOC   = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.OTHER;
x = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA;
% x is previous size, y is current size.
if isempty(x)
  ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA = y;
else
  ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1) = max([x(1),y(1)]);
  ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(2) = max([x(2),y(2)]);
end
return
end

function y = ForZeros(Inputs,callerstr)
global ADIGATOR ADIGATORFORDATA
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
if isa(Inputs{1},'cada')
  ADIGATOR.VARINFO.LASTOCC(Inputs{1}.id,1) =...
    ADIGATOR.VARINFO.COUNT;
end
if length(Inputs) ==2 && isa(Inputs{2},'cada')
  ADIGATOR.VARINFO.LASTOCC(Inputs{2}.id,1) =...
    ADIGATOR.VARINFO.COUNT;
end

% Get the OverMapped Size.
INNERLOC   = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.OTHER;
FMrow = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1);
FNcol = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(2);


% --------------------Build Function Properties-----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if isinf(FMrow)
  if length(Inputs) == 1
    vecDim = [Inputs{1}.func.name,'(1)'];
  else
    vecDim = Inputs{1}.func.name;
  end
  FMrow = 1; yvec = 1;
elseif isinf(FNcol)
  if length(Inputs) == 1
    vecDim = [Inputs{1}.func.name,'(2)'];
  else
    vecDim = Inputs{2}.func.name;
  end
  FNcol = 1; yvec = 2;
else
  yvec = 0;
end
switch callerstr
  case 'zeros'
    y.func.value = zeros(FMrow,FNcol);
  case 'ones'
    y.func.value = ones(FMrow,FNcol);
  case 'eye'
    y.func.value = eye(FMrow,FNcol);
  case 'true'
    y.func.value = true(FMrow,FNcol);
  case 'false'
    y.func.value = false(FMrow,FNcol);
end

if yvec == 1
  fprintf(fid,[indent,funcstr,' = ',callerstr,'(',vecDim,',%1.0f);\n'],FNcol);
elseif yvec == 2
  fprintf(fid,[indent,funcstr,' = ',callerstr,'(%1.0f,',vecDim,');\n'],FMrow);  
else
  fprintf(fid,[indent,funcstr,' = ',callerstr,'(%1.0f,%1.0f);\n'],FMrow,FNcol);
end

ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;

y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
return
end