function y = length(x)
% CADA overloaded LENGTH
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
NUMvod = ADIGATOR.NVAROFDIFF;
indent = ADIGATOR.PRINT.INDENT;

if ADIGATOR.FORINFO.FLAG
  IncreaseForSizeCount();
  if ADIGATOR.RUNFLAG == 2
    y = ForLength(x);
    return
  end
end
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end

y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func  = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',...
  max(x.func.size));
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

if ADIGATOR.PRINT.FLAG
  fprintf(fid,[indent,funcstr,' = length(',x.func.name,');\n']);
end
ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;

ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.RUNFLAG == 1 && ADIGATOR.FORINFO.FLAG
  AssignForSizes(y.func.value);
end
return
end

function IncreaseForSizeCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.SIZE =...
  ADIGATORFORDATA(INNERLOC).COUNT.SIZE +1;
return
end

function AssignForSizes(xSize)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Scount    = ADIGATORFORDATA(INNERLOC).COUNT.SIZE;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
xSize(xSize ==0) = NaN;
% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES = xSize;
else
  ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(1,ITERCOUNT) = xSize;
end

return
end

function y = ForLength(x)
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
Scount   = ADIGATORFORDATA(INNERLOC).COUNT.SIZE;

NUMvod    = ADIGATOR.NVAROFDIFF;
indent    = ADIGATOR.PRINT.INDENT;
fid       = ADIGATOR.PRINT.FID;
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;

y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func  = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

if iscell(ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES)
  IndName = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES{1};
  DepFlag = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES{3}(1);
  IndFlag = 1;
else
  IndFlag = 0;
  outSize = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES;
end

if ~IndFlag
  if isinf(outSize)
    fprintf(fid,[indent,funcstr,' = length(',x.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = %1.0f;\n'],outSize);
  end
elseif DepFlag
  fprintf(fid,[indent,funcstr,' = ',IndName,'(',CountName,');\n']);
else
  fprintf(fid,[indent,funcstr,' = ',IndName,';\n']);
end

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
return
end