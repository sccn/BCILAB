function y = colon(varargin)
% CADA overloaded version of function COLON.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
if ADIGATOR.FORINFO.FLAG > 0
  IncreaseForOtherCount();
end

if nargin == 2
  % y = J:K
  if ADIGATOR.EMPTYFLAG
    y = cadaEmptyEval(varargin{2},varargin{2});
    return
  elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    y = ForColon(varargin{1},varargin{2});
    return
  end
  if isa(varargin{1},'cada')
    if ~isempty(varargin{1}.func.value)
      J = varargin{1}.func.value;
      Jstr = varargin{1}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
    else
      error('colon my not be used with strictly symbolic CADA objects')
    end
  else
    J = varargin{1}(1);
    Jstr = sprintf('%1.0f',J);
  end
  if isa(varargin{2},'cada')
    if ~isempty(varargin{2}.func.value)
      K = varargin{2}.func.value;
      Kstr = varargin{2}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
    else
      error('colon my not be used with strictly symbolic CADA objects')
    end
  else
    K = varargin{2}(1);
    Kstr = sprintf('%1.0f',K);
  end
  
  if isinf(K) && J==1
    vec = 1;
  elseif isinf(J) || isinf(K)
    error(['If N is vectorized dimension, may only use colon as 1:N, ',...
      '1:1:N, or N:-1:1']);
  else
    vec = 0;
    I = J:K;
  end
  Istr = [Jstr,':',Kstr];
  
elseif nargin == 3
  % y = J:D:K
  if ADIGATOR.EMPTYFLAG
    y = cadaEmptyEval(varargin{1},varargin{2},varargin{3});
    return
  elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    y = ForColon(varargin{1},varargin{2},varargin{3});
    return
  end
  if isa(varargin{1},'cada')
    if ~isempty(varargin{1}.func.value)
      J = varargin{1}.func.value;
      Jstr = varargin{1}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
    else
      error('colon my not be used with strictly symbolic CADA objects')
    end
  else
    J = varargin{1}(1);
    Jstr = sprintf('%1.0f',J);
  end
  if isa(varargin{2},'cada')
    if ~isempty(varargin{2}.func.value)
      D = varargin{2}.func.value;
      Dstr = varargin{2}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
    else
      error('colon my not be used with strictly symbolic CADA objects')
    end
  else
    D = varargin{2}(1);
    Dstr = sprintf('%1.0f',D);
  end  
  if isa(varargin{3},'cada')
    if ~isempty(varargin{3}.func.value)
      K = varargin{3}.func.value;
      Kstr = varargin{3}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{3}.id,1) = ADIGATOR.VARINFO.COUNT;
    else
      error('colon my not be used with strictly symbolic CADA objects')
    end
  else
    K = varargin{3}(1);
    Kstr = sprintf('%1.0f',K);
  end  
  if (isinf(K) && J==1 && D == 1) 
    vec=1;
  elseif (isinf(J) && D==-1 && K==1)
    vec=-1;
  elseif isinf(J) || isinf(D) || isinf(K)
    error(['If N is vectorized dimension, may only use colon as 1:N, ',...
      '1:1:N, or N:-1:1']);
  else
    vec=0;
    I = J:D:K;
  end
  Istr = [Jstr,':',Dstr,':',Kstr];
else
  error('Not enough input arguments.')
end

y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
if vec == 1
  y.func = struct('name',funcstr,'size',[1 Inf],'zerolocs',[],'value',[],...
    'OnetoN',1);
elseif vec
  y.func = struct('name',funcstr,'size',[1 Inf],'zerolocs',[],'value',[]);
else
  y.func = struct('name',funcstr,'size',size(I),'zerolocs',[],'value',I);
end

y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

if PFLAG == 1
  fprintf(fid,[indent,funcstr,' = ',Istr,';\n']);
end

ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 1
  ForOtherData(y);
end

return
end

function IncreaseForOtherCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.OTHER = ADIGATORFORDATA(INNERLOC).COUNT.OTHER+1;
return
end

function ForOtherData(y)
global ADIGATOR ADIGATORFORDATA
% All we need to look for here is if the size changes, will store the 
%overmapped size
yNcol      = y.func.size(2);
INNERLOC   = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.OTHER;

if isempty(ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA)
  ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA = yNcol;
else
  xNcol = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA;
  if xNcol ~= yNcol
    ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA = max([xNcol,yNcol]);
  end
end

return
end

function y = ForColon(varargin)
global ADIGATOR ADIGATORFORDATA
NUMvod      = ADIGATOR.NVAROFDIFF;
fid         = ADIGATOR.PRINT.FID;
indent      = ADIGATOR.PRINT.INDENT;
NDstr       = sprintf('%1.0f',ADIGATOR.DERNUMBER);
INNERLOC    = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.OTHER;
yNcol       = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1);


y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func = struct('name',funcstr,'size',[1 yNcol],'zerolocs',[],'value',[]);



y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
tempfuncstr = ['cada',NDstr,'tempf1'];
if nargin == 2
  % y = J:K
  jNum = []; kNum = [];
  if isa(varargin{1},'cada')
      Jstr = varargin{1}.func.name;
      ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
      if ~isempty(varargin{1}.func.value)
        jNum = varargin{1}.func.value(1);
      end
  else
    jNum = varargin{1}(1);
    Jstr = sprintf('%1.0f',jNum);
  end
  if isa(varargin{2},'cada')
    Kstr = varargin{2}.func.name;
    ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
    if ~isempty(varargin{2}.func.value)
      kNum = varargin{2}.func.value(1);
    end
  else
    kNum = varargin{2}(1);
    Kstr = sprintf('%1.0f',kNum);
  end
  Istr = [Jstr,':',Kstr];
  if ~isempty(jNum) && ~isempty(kNum)
    fprintf(fid,[indent,funcstr,' = ',Istr,';\n']);
    if ~isinf(kNum)
      y.func.value = jNum:kNum;
    end
  else
    fprintf(fid,[indent,tempfuncstr,' = zeros(1,%1.0f);\n'],y.func.size(2));
    if strcmp(Jstr,'1')
      fprintf(fid,[indent,tempfuncstr,'(1:',Kstr,') = ',Istr,';\n']);
    else
      fprintf(fid,[indent,tempfuncstr,'(1:',Kstr,'-',Jstr,'+1) = ',Istr,';\n']);
    end
    fprintf(fid,[indent,funcstr,' = ',tempfuncstr,';\n']);
  end
  
elseif nargin == 3
  % y = J:D:K
  jNum = []; dNum = []; kNum = [];
  if isa(varargin{1},'cada')
    Jstr = varargin{1}.func.name;
    ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
    if ~isempty(varargin{1}.func.value)
      jNum = varargin{1}.func.value(1);
    end
  else
    jNum = varargin{1}(1);
    Jstr = sprintf('%1.0f',jNum);
  end
  if isa(varargin{2},'cada')
    Dstr = varargin{2}.func.name;
    ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
    if ~isempty(varargin{2}.func.value)
      dNum = varargin{2}.func.value(1);
    end
  else
    dNum = varargin{2}(1);
    Dstr = sprintf('%1.0f',dNum);
  end
  if isa(varargin{3},'cada')
    Kstr = varargin{3}.func.name;
    ADIGATOR.VARINFO.LASTOCC(varargin{3}.id,1) = ADIGATOR.VARINFO.COUNT;
    if ~isempty(varargin{3}.func.value)
      kNum = varargin{3}.func.value(1);
    end
  else
    kNum = varargin{3}(1);
    Kstr = sprintf('%1.0f',kNum);
  end
  Istr = [Jstr,':',Dstr,':',Kstr];
  if ~isempty(jNum) && ~isempty(dNum) && ~isempty(kNum)
    fprintf(fid,[indent,funcstr,' = ',Istr,';\n']);
    if ~isinf(jNum) && ~isinf(kNum)
      y.func.value = jNum:dNum:kNum;
    end
  else
    fprintf(fid,[indent,tempfuncstr,' = zeros(1,%1.0f);\n'],y.func.size(2));
    fprintf(fid,[indent,tempfuncstr,'(1:(',Kstr,'-',Jstr,'+1)/',Dstr,') = ',Istr,';\n']);
    fprintf(fid,[indent,funcstr,' = ',tempfuncstr,';\n']);
  end
else
  error('Not enough input arguments.')
end

ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT           = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
return
end