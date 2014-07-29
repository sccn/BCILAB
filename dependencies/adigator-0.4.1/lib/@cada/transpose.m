function y = transpose(x)
% CADA overloaded version of function TRANSPOSE.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;

if ADIGATOR.FORINFO.FLAG
  IncreaseForTransposeCount();
  if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    [y,flag,x] = ForTranspose(x);
    if flag; return; end
  end
end
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end
xMrow = x.func.size(1); xNcol = x.func.size(2);
FMrow = xNcol; FNcol = xMrow;
if isinf(FMrow); yvec = 1; elseif isinf(FNcol); yvec = 2; else yvec = 0; end

% --------------------Build Function Properties-----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func  = struct('name',funcstr,'size',[FMrow FNcol],...
  'zerolocs',x.func.zerolocs,'value',x.func.value.');

if isfield(x.func,'logicref')
  y.func.logicref = x.func.logicref([2 1]);
end
if yvec == 1; FMrow = 1; elseif yvec == 2; FNcol = 1; end

% ------------------Build Derivative Properties-----------------------%
xref = zeros(FNcol,FMrow); xref(:) = 1:FMrow*FNcol;
yref = xref.';
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    if FMrow == 1 || FNcol == 1
      % Linear Derivative Index does not change.
      y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
      if DPFLAG == 1
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
    else
      % Linear Derivative Index changes, need to print out the mapping.
      nv    = ADIGATOR.VAROFDIFF(Vcount).usize;
      xrows = x.deriv(Vcount).nzlocs(:,1); 
      xcols = x.deriv(Vcount).nzlocs(:,2);
      dx = sparse(xrows,xcols,1:length(xrows),xMrow*xNcol,nv);
      [yrows,ycols,yind] = find(dx(yref(:),:));
      if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
      y.deriv(Vcount).nzlocs = [yrows, ycols];
      if DPFLAG == 1
        TDind1 = cadaindprint(yind);
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',TDind1,');\n']);
      end
    end
    
  end
end

% -----------------------Function Printing ---------------------------%
if PFLAG == 1
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'.'';\n']);
end
ADIGATOR.VARINFO.LASTOCC([y.id x.id],1)  = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                   = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG
  AssignForTransposeData(y,x);
end
return
end

function IncreaseForTransposeCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.TRANSPOSE = ...
  ADIGATORFORDATA(INNERLOC).COUNT.TRANSPOSE + 1;
return
end

function AssignForTransposeData(y,x)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Tcount    = ADIGATORFORDATA(INNERLOC).COUNT.TRANSPOSE;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;

x.func.size(isinf(x.func.size)) = 1;
y.func.size(isinf(y.func.size)) = 1;
% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES =...
    y.func.size.';
else
  ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES(:,ITERCOUNT) = ...
    y.func.size.';
end

% Variable OverMapping
if ~isa(x,'cada'); x = class(x,'cada'); end
if isempty(ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS)
  % First Call
  ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{1} = y;
  if ~isempty(x.id) && ADIGATOR.VARINFO.NAMELOCS(x.id)
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2} = x;
    end
  else
    ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2} = x;
  end
else
  yOver = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{1};
  ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{1} = cadaUnionVars(y,yOver);
  xOver = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2};
  if isa(xOver,'cada')
    ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2} = cadaUnionVars(x,xOver);
  end
end
return
end

function [y,flag,x] = ForTranspose(x)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Tcount    = ADIGATORFORDATA(INNERLOC).COUNT.TRANSPOSE;

xOver = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{2};
% Check that X is overmapped properly
x = cadaPrintReMap(x,xOver,x.id);

if isempty(ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES)
  % If sizes doesnt change, we dont need to do this.
  flag = 0; y = []; return
else
  flag = 1;
end
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NUMvod    = ADIGATOR.NVAROFDIFF;
NDstr     = sprintf('%1.0d',ADIGATOR.DERNUMBER);
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;

% Get OverMapped Variables
yOver = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).VARS{1};


% Build Y
y = yOver;
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func.name = funcstr;


% ----------------Print Out Derivatives---------------------------------- %
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      
      IndName = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).INDICES{Vcount,1};
      DepFlag = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).INDICES{Vcount,3}(1);
      SpFlag  = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).INDICES{Vcount,3}(2);
      if DepFlag
        % Indices are dependent upon this loop
        IndRef = [IndName,'(:,',CountName,')'];
      else
        % Indices are independent of this loop
        IndRef = IndName;
      end
      if SpFlag
        % Indices are sparse
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
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
% RowInds are Y's Row Sizes, ColInds are Y's Col Sizes
RowInds = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{1,1};
ColInds = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{2,1};
if ~isempty(RowInds) && ~isempty(ColInds)
  RowDepFlag = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{1,3}(1);
  ColDepFlag = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{2,3}(1);
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
  yStr = [funcstr,'(1:',RowRef,',1:',ColRef,')'];
  xStr = [x.func.name,'(1:',ColRef,',1:',RowRef,')'];
elseif ~isempty(RowInds)
  RowDepFlag = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{1,3}(1);
  if RowDepFlag
    RowRef = [RowInds,'(',CountName,')'];
  else
    RowRef = RowInds;
  end
  yStr = [funcstr,'(1:',RowRef,',:)'];
  xStr = [x.func.name,'(:,1:',RowRef,')'];
elseif ~isempty(ColInds)
  ColDepFlag = ADIGATORFORDATA(INNERLOC).TRANSPOSE(Tcount).SIZES{2,3}(1);
  if ColDepFlag
    ColRef = [ColInds,'(',CountName,')'];
  else
    ColRef = ColInds;
  end
  yStr = [funcstr,'(:,1:',ColRef,')'];
  xStr = [x.func.name,'(1:',ColRef,',:)'];
end
% Print out the Temp TRANSPOSE
fprintf(fid,[indent,TF1,' = ',xStr,'.'';\n']);
% Print out Y
fprintf(fid,[indent,funcstr,' = zeros(%1.0d,%1.0d);\n'],y.func.size(1),y.func.size(2));
fprintf(fid,[indent,yStr,' = ',TF1,';\n']);

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
  y = class(y,'cada');
end
return
end
