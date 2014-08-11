function y = vertcat(varargin)
% CADA overloaded version of VERTCAT
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);

if ADIGATOR.FORINFO.FLAG
  IncreaseForVertcatCount();
end
if ADIGATOR.EMPTYFLAG
  % Need to send this to cadaEmptyEval, build the string to call it.
  InputStr = cell(1,nargin);
  for Icount = 1:nargin
    InputStr{Icount} = sprintf('varargin{%1.0d},',Icount);
  end
  InputStr = cell2mat(InputStr);
  ForOpCall = ['cadaEmptyEval(',InputStr(1:end-1),');'];
  y = eval(ForOpCall);
  return
elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
  [y,flag,varargin] = ForVertcat(varargin);
  % flag = 0 => don't need to do anything special for this horzcat
  if flag; return; end
end

% ------------Parse Input Sizes/Make Numerics to Overloaded-------------- %
if PFLAG; InputStr = cell(1,nargin); end
AllNumericFlag = 1;
iMrows = zeros(1,nargin); RowSum = zeros(1,nargin);
for Icount = 1:nargin
  x = varargin{Icount};
  if isa(x,'cada')
    xMrow = x.func.size(1); xNcol = x.func.size(2);
    if isempty(x.func.value) && xMrow > 0 && xNcol > 0
      AllNumericFlag = 0;
    end
  elseif isnumeric(x)
    x = Num2Overloaded(x);
    xMrow = x.func.size(1); xNcol = x.func.size(2);
    varargin{Icount} = x;
  else
    error('??? Invalid input to VERTCAT')
  end
  if Icount == 1; 
    yNcol = xNcol;
  elseif yNcol ~= xNcol && yNcol > 0 && xNcol > 0
    error('CAT arguments dimensions are not consistent.')
  elseif yNcol == 0
    yNcol = xNcol;
  end
  iMrows(Icount) = xMrow;
  RowSum(Icount) = sum(iMrows);
  if PFLAG; InputStr{Icount} = [x.func.name,';']; end
end

% -----------------Begin Parsing Function Information-------------------- %
y.id    = ADIGATOR.VARINFO.COUNT;
yMrow   = RowSum(nargin);
[funcstr,DPFLAG] = cadafuncname();
y.func  = struct('name',funcstr,'size',[yMrow yNcol],'zerolocs',[],...
  'value',[]);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if isinf(yNcol)
  yvec = 1; yNcol = 1;
else
  yvec = 0;
end
if isinf(yMrow)
  error('Cannot catonate in vectorized dimension')
end
if AllNumericFlag
  % All inputs have a known numeric value - will be no derivatives
  y.func.value = zeros(yMrow,yNcol);
  % ---Make Another Sweep through Inputs to get Numeric y values---
  y.func.value(1:iMrows(1),:) = varargin{1}.func.value;
  for Icount = 2:nargin
    y.func.value(RowSum(Icount-1)+1:RowSum(Icount),:) =...
      varargin{Icount}.func.value;
  end
else
  % At least one of the inputs is symbolic
  % Build a Reference Matrix for each Input
  % y indices corresponding to input i can be found by iRefs{i}(i_indices)
  yRef = zeros(yMrow,yNcol); yRef(:) = 1:yMrow*yNcol;
  iRefs = cell(1,nargin); iRefs{1} = yRef(1:iMrows(1),:);
  for Icount = 2:nargin
    iRefs{Icount} = yRef(RowSum(Icount-1)+1:RowSum(Icount),:);
  end
  % ---Make Another Sweep through inputs to gain function sparsity and
  % number of derivatives---
  yTemp       = true(yMrow,yNcol);
  iNumDerivs  = zeros(NUMvod,nargin); 
  iDerivSum   = zeros(NUMvod,nargin);

  for Icount = 1:nargin
    x = varargin{Icount};
    % Function Sparsity
    yTemp(iRefs{Icount}(x.func.zerolocs)) = false;
    % Number of Derivatives
    for Vcount = 1:NUMvod
      if ~isempty(x.deriv(Vcount).nzlocs)
        iNumDerivs(Vcount,Icount) = size(x.deriv(Vcount).nzlocs,1);
      end
      iDerivSum(Vcount,Icount) = sum(iNumDerivs(Vcount,:));
    end
  end
  
  if nnz(yTemp) < numel(yTemp)
    % function sparsity
    y.func.zerolocs = find(~yTemp(:));
  end
  
  % -------------------------Build Derivative---------------------------- %
  for Vcount = 1:NUMvod
    yNumDerivs = iDerivSum(Vcount,nargin);
    if yNumDerivs
      % At least one input has derivatives wrt this variable
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      dyMat = zeros(yNumDerivs,3);
      % First Pass through to get DY
      for Icount = 1:nargin
        x = varargin{Icount};
        if ~isempty(x.deriv(Vcount).nzlocs)
          xrows = x.deriv(Vcount).nzlocs(:,1);
          xcols = x.deriv(Vcount).nzlocs(:,2);

          xrows = iRefs{Icount}(xrows);
          if size(xrows,2) > 1; xrows = xrows.'; end
          if Icount == 1
            dyMat(1:iDerivSum(Vcount,Icount),:) = [xrows xcols ones(size(xrows))];
          else
            dyMat(iDerivSum(Vcount,Icount-1)+1:...
              iDerivSum(Vcount,Icount),:) = [xrows xcols Icount*ones(size(xrows))];
          end
        end
      end
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      dy = sparse(dyMat(:,1),dyMat(:,2),dyMat(:,3),yMrow*yNcol,nv);

      [yrows ycols] = find(dy);
      if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
      y.deriv(Vcount).nzlocs = [yrows ycols];

      if DPFLAG
        % Need to make another sweep to get the mapping of all of the
        % derivatives and print them out.
        TD1 = ['cada',NDstr,'td1'];
        if yvec
          fprintf(fid,[indent,TD1,' = zeros(size(',varargin{1}.func.name,',2),%1.0d);\n'],yNumDerivs);
        else
          fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],yNumDerivs);
        end
        dy = sparse(yrows,ycols,1:yNumDerivs,yMrow*yNcol,nv);
        for Icount = 1:nargin
          x = varargin{Icount};
          if ~isempty(x.deriv(Vcount).nzlocs)
            xinds  = nonzeros(dy(iRefs{Icount}(:),:));
            TDind1 = cadaindprint(xinds);
            if yvec
              fprintf(fid,[indent,TD1,'(:,',TDind1,') = ',...
                x.deriv(Vcount).name,';\n']);
            else
              fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
                x.deriv(Vcount).name,';\n']);
            end
          end
        end
        fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
      end
    end
  end
end

if PFLAG
  % Print out the Function
  InputStr = cell2mat(InputStr);
  fprintf(fid,[indent,funcstr,' = [',InputStr(1:end-1),'];\n']);
end

ADIGATOR.VARINFO.LASTOCC(y.id,1)   = ADIGATOR.VARINFO.COUNT;
for Icount = 1:nargin
  x = varargin{Icount};
  ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
end
y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ADIGATOR.FORINFO.FLAG
  AssignForVertcatData(yNcol,iMrows,y,varargin);
end
return
end

function IncreaseForVertcatCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.VERTCAT =...
  ADIGATORFORDATA(INNERLOC).COUNT.VERTCAT + 1;
return
end

function AssignForVertcatData(yNcol,iMrows,y,Inputs)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
VEcount   = ADIGATORFORDATA(INNERLOC).COUNT.VERTCAT;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
NUMinput  = length(iMrows);

% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).SIZES = [yNcol,iMrows].';
else
  ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).SIZES(1:1+NUMinput,ITERCOUNT)...
    = [yNcol,iMrows].';
end

% Variable OverMapping
if isempty(ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS)
  % First Call
  ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS = cell(1+NUMinput,1);
  ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{1} = y;
  for Icount = 1:NUMinput
    x = Inputs{Icount};
    if ~isempty(x.id) && ADIGATOR.VARINFO.NAMELOCS(x.id)
      OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
      StartCount = ADIGATORFORDATA(OUTERLOC).START;
      EndCount   = ADIGATORFORDATA(OUTERLOC).END;
      xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
      xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
      if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
        ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} = xOverLoc1;
      elseif xOverLoc2 && ...
          any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
        ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} = xOverLoc2;
      else
        ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} = x;
      end
    else
      ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} = x;
    end
  end
elseif length(ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS) ~= NUMinput+1
  % Can't think of a case for this ever happening, but who knows
  error('??? Number of inputs to VERTCAT cannot change size in a for loop.')
else
  % Output Variable - y
  yOver = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{1};
  ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{1} = cadaUnionVars(yOver,y);
  for Icount = 1:NUMinput
    % Input Variables
    x     = Inputs{Icount};
    xOver = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1};
    if isa(xOver,'cada')
      ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} =...
        cadaUnionVars(xOver,x);
    else
      xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
      if xOverLoc1 && xOver ~= xOverLoc1
        ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1} = xOverLoc1;
      end
    end
  end
end
return
end

function [y,OutFlag,Inputs] = ForVertcat(Inputs)
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
VEcount  = ADIGATORFORDATA(INNERLOC).COUNT.VERTCAT;
if isempty(ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).SIZES)
  % Sizes dont change - dont need to do anything special
  TVcount = 0;
  NUMinputs = length(Inputs);
  for Icount = 1:NUMinputs
    iOver = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1};
    % Check for numeric input
    if isnumeric(Inputs{Icount})
      TVcount = TVcount+1;
      Inputs{Icount} = Num2Overloaded(Inputs{Icount},TVcount);
    end
    % Check to make sure input and Overmapped inputs match
    Inputs{Icount} = cadaPrintReMap(Inputs{Icount},iOver,Inputs{Icount}.id);
  end
  OutFlag = 0; y = []; return
else
  OutFlag = 1;
end
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0d',ADIGATOR.DERNUMBER);
NUMvod    = ADIGATOR.NVAROFDIFF;
NUMinputs = length(Inputs);

% ---Get OverMapped Inputs--- %
yOver   = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{1};
TVcount = 0;
for Icount = 1:NUMinputs
  iOver = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).VARS{Icount+1};
  % Check for numeric input
  if isnumeric(Inputs{Icount})
    TVcount = TVcount+1;
    Inputs{Icount} = Num2Overloaded(Inputs{Icount},TVcount);
  end
  % Check to make sure input and Overmapped inputs match
  Inputs{Icount} = cadaPrintReMap(Inputs{Icount},iOver,Inputs{Icount}.id);
end

y                 = yOver;
y.id              = ADIGATOR.VARINFO.COUNT;
[funcstr DPFLAG]  = cadafuncname();
y.func.name       = funcstr;
if isinf(y.func.size(2)); yvec = 1; else yvec = 0; end


% -------------------Print Out Derivatives------------------------------- %
for Vcount = 1:NUMvod
  if ~isempty(y.deriv(Vcount).nzlocs)
    y.deriv(Vcount).name = cadadername(funcstr,Vcount);
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      if yvec
        fprintf(fid,[indent,TD1,' = zeros(size(',Inputs{1}.func.name,',2),%1.0d);\n'],size(y.deriv(Vcount).nzlocs,1));
      else
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
      end
      for Icount = 1:NUMinputs
        if ~isempty(Inputs{Icount}.deriv(Vcount).nzlocs)
          IndName = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).INDICES{Vcount,Icount,1};
          DepFlag = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).INDICES{Vcount,Icount,3}(1);
          SpFlag  = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).INDICES{Vcount,Icount,3}(2);
          if DepFlag % Dependent on this loop
            IndRef = [IndName,'(:,',CountName,')'];
          else % Independent of this loop
            IndRef = IndName;
          end
          if yvec
            if SpFlag % Sparse References
              fprintf(fid,[indent,TD1,'(:,logical(',IndRef,')) = ',...
                Inputs{Icount}.deriv(Vcount).name,'(:,nonzeros(',IndRef,'));\n']);
            else % Dense References
              fprintf(fid,[indent,TD1,'(:,',IndRef,') = ',...
                Inputs{Icount}.deriv(Vcount).name,'(:,',IndRef,');\n']);
            end
          else
            if SpFlag % Sparse References
              fprintf(fid,[indent,TD1,'(logical(',IndRef,')) = ',...
                Inputs{Icount}.deriv(Vcount).name,'(nonzeros(',IndRef,'));\n']);
            else % Dense References
              fprintf(fid,[indent,TD1,'(',IndRef,') = ',...
                Inputs{Icount}.deriv(Vcount).name,'(',IndRef,');\n']);
            end
          end
        end
      end
      fprintf(fid,[indent,y.deriv(Vcount).name,' = ',TD1,';\n']);
    end
  end
end

% ---------------------Print Out Function-------------------------------- %
InputString = cell(1,NUMinputs);
for Icount = 1:NUMinputs
  IndName = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).SIZES{Icount,1};
  if ~isempty(IndName)
    DepFlag = ADIGATORFORDATA(INNERLOC).VERTCAT(VEcount).SIZES{Icount,3}(1);
    if DepFlag
      InputString{Icount} = [Inputs{Icount}.func.name,'(1:',IndName,'(',CountName,'),:);'];
    else
      InputString{Icount} = [Inputs{Icount}.func.name,'(1:',IndName,',:);'];
    end
  else
    InputString{Icount} = [Inputs{Icount}.func.name,';'];
  end
end

TF1 = ['cada',NDstr,'tempf1'];
InputString = cell2mat(InputString);
fprintf(fid,[indent,TF1,' = [',InputString(1:end-1),'];\n']);
if yvec
  fprintf(fid,[indent,funcstr,' = zeros(%1.0d,size(',Inputs{1}.func.name,',2));\n'],y.func.size(1));
else
  fprintf(fid,[indent,funcstr,' = zeros(%1.0d,%1.0d);\n'],y.func.size(1),y.func.size(2));
end
fprintf(fid,[indent,funcstr,'(1:size(',TF1,',1),:) = ',TF1,';\n']);

ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;
for Icount = 1:NUMinputs
  x = Inputs{Icount};
  ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
end
if ~isa(y,'cada')
  y = class(y,'cada');
end
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
return
end

function y = Num2Overloaded(x)
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
if ADIGATOR.PRINT.FLAG
  yname   = cadamatprint(x);
else
  yname = [];
end
y.id = [];
y.func  = struct('name',yname,'size',size(x),'zerolocs',[],'value',x);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

if nnz(x) < numel(x)
  [xrows,xcols] = find(x);
  if size(xrows,2) > 1; xrows = xrows.'; xcols = xcols.'; end
  y.func.zerolocs = [xrows,xcols];
end
y = class(y,'cada');

return
end