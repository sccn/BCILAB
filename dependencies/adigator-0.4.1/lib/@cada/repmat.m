function y = repmat(x,varargin)
% CADA overloaded version of REPMAT
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;

if ADIGATOR.FORINFO.FLAG
  IncreaseForRepmatCount();
  if ADIGATOR.RUNFLAG == 2
    [y,flag,x] = ForRepmat(x,varargin);
    if flag; return; end
  end
end
if ADIGATOR.EMPTYFLAG
  % Need to send to cadaEmptyEval, build the string to call it
  InputStr = cell(1,nargin);
  InputStr{1} = 'x,';
  for Icount = 2:nargin
    InputStr{Icount} = sprintf('varargin{%1.0f},',Icount-1);
  end
  InputStr = cell2mat(InputStr);
  ForOpCall = ['cadaEmptyEval(',InputStr(1:end-1),');'];
  y = eval(ForOpCall);
  return
end
% ------------------Parse Sizing Inputs---------------------------------- %
if nargin == 1
  error('Requires at least 2 inputs.')
elseif nargin == 2
  if isa(varargin{1},'cada')
    if varargin{1}.func.size(1)*varargin{1}.func.size(2) == 1
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value;
        repNcol = repMrow;
      else
        error('Cannot REPMAT a purely symbolic dimension.')
      end
    elseif varargin{1}.func.size(1)*varargin{1}.func.size(2) == 2
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value(1);
        repNcol = varargin{1}.func.value(2);
        if isinf(repMrow); 
          vecDimStr = [varargin{1}.func.name,'(1)'];
        elseif isinf(repNcol); 
          vecDimStr = [varargin{1}.func.name,'(2)']; 
        end
      else
        error('Cannot REPMAT a purely symbolic dimension.')
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{1})
    if length(varargin{1}) == 1
      repMrow = varargin{1};
      repNcol = repMrow;
    elseif length(varargin{1}) == 2
      repMrow = varargin{1};
      repNcol = varargin{2};
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
elseif nargin == 3
  if isa(varargin{1},'cada')
    if varargin{1}.func.size(1)*varargin{1}.func.size(2) == 1
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value;
        if isinf(repMrow)
          vecDimStr = varargin{1}.func.name;
        end
      else
        error('Cannot REPMAT a purely symbolic dimension.')
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{1})
    if length(varargin{1}) == 1
      repMrow = varargin{1};
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
  if isa(varargin{2},'cada')
    if varargin{2}.func.size(1)*varargin{2}.func.size(2) == 1
      if ~isempty(varargin{2}.func.value)
        repNcol = varargin{2}.func.value;
        if isinf(repNcol)
          vecDimStr = varargin{2}.func.name;
        end
      else
        error('Cannot REPMAT a purely symbolic dimension.')
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{2})
    if length(varargin{2}) == 1
      repNcol = varargin{2};
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
else
  error('Too many input arguments.')
end
% --------------------------Parse X-------------------------------------- %
if isa(x,'cada')
  xMrow = x.func.size(1);
  xNcol = x.func.size(2);
elseif isnumeric(x)
  [xMrow,xNcol] = size(x);
  xtemp.id = [];
  xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
    'value',x);
  if PFLAG
    xtemp.func.name = cadamatprint(x);
  end
  xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  x = xtemp;
else
  y = repmat(x,repMrow,repNcol);
  return
end

% --------------------Build Function Properties-----------------------%
FMrow = xMrow*repMrow;
FNcol = xNcol*repNcol;
y.id             = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);
if isinf(FMrow) && isinf(xMrow) && repMrow == 1 && ~isinf(FNcol)
  xMrow = 1; rvec = 0; 
elseif isinf(FMrow) && isinf(repMrow) && xMrow == 1 && ~isinf(FNcol)
  repMrow = 1; rvec = 1;
elseif isinf(FNcol) && isinf(xNcol) && repNcol == 1 && ~isinf(FMrow)
  xNcol = 1; rvec = 0;
elseif isinf(FNcol) && isinf(repNcol) && xNcol == 1 && ~isinf(FMrow)
  xNcol = 1; rvec = 2;
elseif ~isinf(FMrow) && ~isinf(FNcol)
  rvec = 0;
else
  error(['Invalid vectorized repmat: if the vectorized dimension is N',...
    ' then the result of the repmat must be of size N by n, or n by N,',...
    ' where n is a non-vectorized known dimension'])
end

findices = zeros(xMrow,xNcol);
findices(:) = 1:xMrow*xNcol;
findices    = repmat(findices,repMrow,repNcol);

% function sparsity
if ~isempty(x.func.zerolocs)
  xtemp = true(xMrow,xNcol); xtemp(x.func.zerolocs) = false;
  ytemp = repmat(xtemp,repMrow,repNcol);
  y.func.zerolocs = find(~ytemp(:));
end
if isfield(x.func,'logical')
  y.func.logical = [];
end
% function numerics (if exist)
y.func.value = repmat(x.func.value,repMrow,repNcol);

% ------------------Build Derivative Properties-----------------------%
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    nzx = size(x.deriv(Vcount).nzlocs,1);
    nv  = ADIGATOR.VAROFDIFF(Vcount).usize;
    dx  = sparse(x.deriv(Vcount).nzlocs(:,1),x.deriv(Vcount).nzlocs(:,2),(1:nzx)',xMrow*xNcol,nv);
    [yrows,ycols,yind] = find(dx(findices(:)',:));
    if size(yrows,2) > 1
      yrows =  yrows';
      ycols = ycols';
    end
    y.deriv(Vcount).nzlocs = [yrows,ycols];
    % -----------Derivative Vectorization/Printing----------------%
    if DPFLAG
      TFind1 = cadaindprint(yind);
      if isinf(FMrow) || isinf(FNcol)
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(:,',TFind1,');\n']);
      else
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',TFind1,');\n']);
      end
    end
  end
end

% -------------Function Vectorization/Printing -----------------------%
if PFLAG
  if isinf(repMrow)
    fprintf(fid,[indent,funcstr,' = repmat(',x.func.name,...
      ',',vecDimStr,',%1.0f);\n'],repNcol);
  elseif isinf(repNcol)
    fprintf(fid,[indent,funcstr,' = repmat(',x.func.name,...
      ',%1.0f,',vecDimStr,');\n'],repMrow);
  else
    fprintf(fid,[indent,funcstr,' = repmat(',x.func.name,...
      ',%1.0f,%1.0f);\n'],repMrow,repNcol);
  end
end

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG
  AssignForRepmatData(y,x);
end
return
end

function IncreaseForRepmatCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.REPMAT = ...
  ADIGATORFORDATA(INNERLOC).COUNT.REPMAT + 1;
return
end

function AssignForRepmatData(y,x)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Rcount    = ADIGATORFORDATA(INNERLOC).COUNT.REPMAT;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
x.func.size(isinf(x.func.size)) = 1;
y.func.size(isinf(y.func.size)) = 1;
% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES =...
    [y.func.size.';x.func.size.'];
else
  ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES(:,ITERCOUNT) = ...
    [y.func.size.';x.func.size.'];
end

% Variable OverMapping
if ~isa(x,'cada'); x = class(x,'cada'); end
if isempty(ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS)
  % First Call
  ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{1} = y;
  if ~isempty(x.id) && ADIGATOR.VARINFO.NAMELOCS(x.id)
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = x;
    end
  else
    ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = x;
  end
else
  yOver = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{1};
  ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{1} = cadaUnionVars(y,yOver);
  xOver = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2};
  if isa(xOver,'cada')
    ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = cadaUnionVars(x,xOver);
  else
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    if xOverLoc1 && xOver ~= xOverLoc1
      ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2} = xOverLoc1;
    end
  end
end
return
end

function [y,flag,x] = ForRepmat(x,Sinputs)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
Rcount    = ADIGATORFORDATA(INNERLOC).COUNT.REPMAT;

xOver = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{2};
% Check that X is overmapped properly
x = cadaPrintReMap(x,xOver,x.id);

if isempty(ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES)
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
      repMstr = Sinputs{1}.func.name;
      RepStr = [repMstr,',',repMstr];
    else
      RepStr = Sinputs{1}.func.name;
    end
    ADIGATOR.VARINFO.LASTOCC(Sinputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{1})
    if length(Sinputs{1}) == 1
      repMrow = Sinputs{1};
      repNcol = repMrow;
    elseif length(Sinputs{1}) == 2
      repMrow = Sinputs{1}(1);
      repNcol = Sinputs{1}(2);
    end
    RepStr = sprintf('%1.0d,%1.0d',repMrow,repNcol);
  end
elseif length(Sinputs) == 2
  if isa(Sinputs{1},'cada')
    repMstr = Sinputs{1}.func.name;
    ADIGATOR.VARINFO.LASTOCC(Sinputs{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{1})
    repMstr = sprintf('%1.0d',Sinputs{1});
  end
  if isa(Sinputs{2},'cada')
    repNstr = Sinputs{2}.func.name;
    ADIGATOR.VARINFO.LASTOCC(Sinputs{2}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(Sinputs{2})
    repNstr = sprintf('%1.0d',Sinputs{2});
  end
  RepStr =[repMstr,',',repNstr];
end

% Get OverMapped Variables
yOver = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).VARS{1};


% Build Y
y = yOver;
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func.name = funcstr;
yMrow = y.func.size(1); xMrow = x.func.size(1);
yNcol = y.func.size(2); xNcol = x.func.size(2);
if isinf(yMrow)
  if isinf(xMrow)
    vecDim = ['size(',x.func.name,',1)'];
  else
    vecDim = repMstr;
  end
  yvec = 1;
elseif isinf(yNcol)
  if isinf(xNcol)
    vecDim = ['size(',x.func.name,',2)'];
  else
    vecDim = repNstr;
  end
  yvec = 1;
else
  yvec = 0;
end


% ----------------Print Out Derivatives---------------------------------- %
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      if yvec
        fprintf(fid,[indent,TD1,' = zeros(',vecDim,',1);\n'],size(y.deriv(Vcount).nzlocs,1));
      else
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
      end
      IndName = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).INDICES{Vcount,1};
      DepFlag = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).INDICES{Vcount,3}(1);
      SpFlag  = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).INDICES{Vcount,3}(2);
      if DepFlag
        % Indices are dependent upon this loop
        IndRef = [IndName,'(:,',CountName,')'];
      else
        % Indices are independent of this loop
        IndRef = IndName;
      end
      if yvec
        if SpFlag
          % Indices are sparse
          fprintf(fid,[indent,TD1,'(:,logical(',IndRef,')) = ',...
            x.deriv(Vcount).name,'(:,nonzeros(',IndRef,'));\n']);
        else
          % Indices are non-sparse
          fprintf(fid,[indent,TD1,' = ',...
            x.deriv(Vcount).name,'(:,',IndRef,');\n']);
        end
      else
        if SpFlag
          % Indices are sparse
          fprintf(fid,[indent,TD1,'(logical(',IndRef,')) = ',...
            x.deriv(Vcount).name,'(nonzeros(',IndRef,'));\n']);
        else
          % Indices are non-sparse
          fprintf(fid,[indent,TD1,' = ',...
            x.deriv(Vcount).name,'(',IndRef,');\n']);
        end
      end
      fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
    end
  end
end

% ----------------------Print Out Function------------------------------- %
TF1 = ['cada',NDstr,'tempf1'];
% Get Our Reference off of X first
RowInds = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{1,1};
ColInds = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{2,1};
if ~isempty(RowInds) && ~isempty(ColInds)
  RowDepFlag = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{1,3}(1);
  ColDepFlag = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{2,3}(1);
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
  RowDepFlag = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{1,3}(1);
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
  ColDepFlag = ADIGATORFORDATA(INNERLOC).REPMAT(Rcount).SIZES{2,3}(1);
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
% Print out the Temp REPMAT
fprintf(fid,[indent,TF1,' = repmat(',xStr,',',RepStr,');\n']);
% Print out Y
if isinf(yMrow)
  fprintf(fid,[indent,funcstr,' = zeros(',vecDim,',%1.0d);\n'],y.func.size(2));
elseif isinf(yNcol)
  fprintf(fid,[indent,funcstr,' = zeros(%1.0d,',vecDim,');\n'],y.func.size(1),y.func.size(2));
else
  fprintf(fid,[indent,funcstr,' = zeros(%1.0d,%1.0d);\n'],y.func.size(1),y.func.size(2));
end
fprintf(fid,[indent,funcstr,'(1:size(',TF1,',1),1:size(',TF1,',2)) = ',TF1,';\n']);

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
  y = class(y,'cada');
end
return
end