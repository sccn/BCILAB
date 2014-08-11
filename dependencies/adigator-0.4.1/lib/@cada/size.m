function varargout = size(x,varargin)
% CADA overloaded SIZE - note: there is a workaround to stop the MATLAB
% workspace from calling this and increasing the variable count which
% amounts to setting the ADIGATOR.OPTIONS.KEYBOARD flag to 1. It is
% advisable to set this to 1 if you are going to mess with this function.
%
% Written by Matthew J. Weinstein
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;

if nargout == 1 && nargin == 1
  % This is the case that Matlab Workspace will call - If you put a
  % keyboard in this section you are going to crash matlab if you have the
  % workspace open.
  Dbstuff = dbstack; CallingFile = Dbstuff(2).file;
  if ADIGATOR.OPTIONS.KEYBOARD
    error(['Cannot use size with nargin = nargout = 1 when you have a',...
      '''keyboard'' within the code - causes an error with the MATLAB WorkSpace.']);
  elseif length(CallingFile) > 16 && strcmp(CallingFile(1:16),'adigatortempfunc')
    % This is being called from user code not a keyboard.
    if ADIGATOR.FORINFO.FLAG
      IncreaseForSizeCount();
      if ADIGATOR.RUNFLAG == 2
        varargout{1} = ForSize(x);
        return
      end
    end
    if ADIGATOR.EMPTYFLAG
      varargout{1} = cadaEmptyEval(x);
      return
    end
    y.id = ADIGATOR.VARINFO.COUNT;
    [funcstr,~] = cadafuncname();
    y.func = struct('name',funcstr,'size',[1 2],'zerolocs',[],'value',x.func.size);
    y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    
    if ADIGATOR.PRINT.FLAG
      fprintf(fid,[indent,funcstr,' = size(',x.func.name,');\n']);
    end
    
    ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    y = class(y,'cada');
    if ADIGATOR.FORINFO.FLAG
      AssignForSizes(y.func.value);
    end
    varargout{1} = y;
  else
    varargout{1} = x.func.size;
  end
elseif nargin == 2 && (nargout == 1 || nargout == 0)
  % Want the size of the first or second dimension
  if ADIGATOR.FORINFO.FLAG
    IncreaseForSizeCount();
    if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
      varargout{1} = ForSize(x,varargin{1});
      return
    end
  end
  if ADIGATOR.EMPTYFLAG
    varargout{1} = cadaEmptyEval(x,varargin{1});
    return
  end
  Dim = varargin{1};
  if isa(Dim,'cada')
    if ~isempty(Dim.func.value)
      ADIGATOR.VARINFO.LASTOCC(Dim.id,1) =...
        ADIGATOR.VARINFO.COUNT;
      DimStr = Dim.func.name;
      Dim = Dim.func.value;
    else
      error('Dimension cannot be purely symbolic.')
    end
  else
    DimStr = sprintf('%1.0d',Dim);
  end
  if numel(Dim) ~= 1 || Dim > 2
    error('Dimension must be either 1 or 2.')
  end
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  if isa(x,'cada')
    ADIGATOR.VARINFO.LASTOCC(x.id,1) =...
      ADIGATOR.VARINFO.COUNT;
    y.func.value = x.func.size(Dim);
    if ADIGATOR.PRINT.FLAG
      % Print it out
      fprintf(fid,[indent,funcstr,' = size(',x.func.name,',',DimStr,');\n']);
    end
  else
    y.func.value = size(x,Dim);
    if ADIGATOR.PRINT.FLAG
      % Print it out
      fprintf(fid,[indent,funcstr,' = %1.0f;\n'],y.func.value);
    end
  end
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  varargout{1} = class(y,'cada');
  if ADIGATOR.FORINFO.FLAG
    AssignForSizes(y.func.value);
  end

elseif nargin == 1 && (nargout == 2 || nargout == 0)
  % 2 Outputs, size of 1st and 2nd dimension
  if ADIGATOR.FORINFO.FLAG
    IncreaseForSizeCount();
    if ADIGATOR.RUNFLAG == 2
      if nargout == 2
        [varargout{1},varargout{2}] = ForSize(x);
      else
        varargout{1} = ForSize(x);
      end
      return
    end
  end
  if ADIGATOR.EMPTYFLAG
    if nargout == 2
      [varargout{1}, varargout{2}] = cadaEmptyEval(x);
    else
      varargout{1} = cadaEmptyEval(x);
    end
    return
  end
  
  y1.id = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.LASTOCC([x.id y1.id],1) =...
    ADIGATOR.VARINFO.COUNT+1;
  [funcstr,~] = cadafuncname();
  y1.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',x.func.size(1));
  y1.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  if nargout == 2
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    y2.id = ADIGATOR.VARINFO.COUNT;
    [funcstr,~] = cadafuncname();
    y2.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',x.func.size(2));
    y2.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    
    if ADIGATOR.PRINT.FLAG
      % Print Out the call
      fprintf(fid,[indent,'[',y1.func.name,', ',y2.func.name,...
        '] = size(',x.func.name,');\n']);
    end
    ADIGATOR.VARINFO.LASTOCC(y2.id,1) = ...
      ADIGATOR.VARINFO.COUNT;
    varargout{2} = class(y2,'cada');
  else
    y1.func.size = [1 2];
    y1.func.value = x.func.size;
  end
  ADIGATOR.VARINFO.LASTOCC(y1.id,1) = ...
    ADIGATOR.VARINFO.COUNT;
  
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  varargout{1} = class(y1,'cada');
  if ADIGATOR.FORINFO.FLAG
    AssignForSizes(x.func.size);
  end
elseif nargout > 2
  error('Too many output arguments')
elseif nargin > 2
  error('Too many input arguments')
end
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

xSize(xSize == 0) = NaN;

% Assign the Sizes
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES = xSize(:);
else
  ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(:,ITERCOUNT) = xSize(:);
end
if length(xSize) == 2 && isinf(xSize(1))
  if any(ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(2,1)~= ...
      ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(2,:))
    error(['Cannot use SIZE with two outputs, where one dimension is ',...
      'vectorized and the other changes within a loop.']);
  end
elseif length(xSize) == 2 && isinf(xSize(2))
  if any(ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(1,1)~= ...
      ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES(1,:))
    error(['Cannot use SIZE outputing two dims, where one dimension is ',...
      'vectorized and the other changes within a loop.']);
  end
end
return
end

function [y, varargout] = ForSize(x,varargin)
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
Scount   = ADIGATORFORDATA(INNERLOC).COUNT.SIZE;

NUMvod    = ADIGATOR.NVAROFDIFF;
indent    = ADIGATOR.PRINT.INDENT;
fid       = ADIGATOR.PRINT.FID;
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;

if isa(x,'cada')
  ADIGATOR.VARINFO.LASTOCC(x.id,1) =...
    ADIGATOR.VARINFO.COUNT;
end
if nargin == 2 && isa(varargin{1},'cada')
  ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) =...
    ADIGATOR.VARINFO.COUNT;
end


% All we are going to do is print out what the sizes are.
if iscell(ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES)
  IndName = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES{1};
  DepFlag = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES{3}(1);
  IndFlag = 1;
else
  IndFlag = 0;
  outSize = ADIGATORFORDATA(INNERLOC).SIZE(Scount).SIZES;
end

if nargin == 2
  % Only one dimension => only one output
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  if ~IndFlag
    if isinf(outSize)
      if isinf(x.func.size(1))
        fprintf(fid,[indent,funcstr,' = size(',x.func.name,',1);\n']);
      else
       fprintf(fid,[indent,funcstr,' = size(',x.func.name,',2);\n']);
      end
    else
      fprintf(fid,[indent,funcstr,' = %1.0f;\n'],outSize);
    end
  elseif DepFlag
    fprintf(fid,[indent,funcstr,' = ',IndName,'(',CountName,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = ',IndName,';\n']);
  end
  ADIGATOR.VARINFO.LASTOCC(y.id,1) =...
    ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');
elseif nargout == 2
  % Two Outputs, both dimensions
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y1.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y1.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);
  y1.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  if ~IndFlag
    if isinf(OutSize(1))
      fprintf(fid,[indent,y.func.name,' = size(',x.func.name,',1);\n']);
    else
      fprintf(fid,[indent,y.func.name,' = %1.0f;\n'],OutSize(1));
    end
    if isinf(OutSize(2))
      fprintf(fid,[indent,y1.func.name,' = size(',x.func.name,',2);\n']);
    else
      fprintf(fid,[indent,y1.func.name,' = %1.0f;\n'],OutSize(2));
    end
  elseif DepFlag
    fprintf(fid,[indent,y.func.name,' = ',IndName,'(1,',CountName,');\n']);
    fprintf(fid,[indent,y1.func.name,' = ',IndName,'(2,',CountName,');\n']);
  else
    fprintf(fid,[indent,y.func.name,' = ',IndName,'(1);\n']);
    fprintf(fid,[indent,y1.func.name,' = ',IndName,'(2);\n']);
  end
  
  ADIGATOR.VARINFO.LASTOCC([y.id y1.id],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');
  varargout{1} = class(y1,'cada');
elseif nargout == 1
  % One Output, both dimensions
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y.func = struct('name',funcstr,'size',[1 2],'zerolocs',[],'value',[]);
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  if ~IndFlag
    if any(isinf(outSize))
      fprintf(fid,[indent,funcstr,' = size(',x.func.name,');\n']);
    else
      fprintf(fid,[indent,funcstr,' = [%1.0f,%1.0f];\n'],outSize(1),outSize(2));
    end
  elseif DepFlag
    fprintf(fid,[indent,funcstr,' = ',IndName,'(:,',CountName,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = ',IndName,';\n']);
  end
  
  
  ADIGATOR.VARINFO.LASTOCC(y.id,1) =...
    ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');
else
  error('??? too many outputs')
end

return
end