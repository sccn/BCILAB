function y = subsref(x,s)
% CADA overloaded version of function SUBSREF.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
ssize = length(s);
for scount = 1:ssize;
  if strcmp(s(scount).type,'()')
    if scount == 1
      % ------------- Overloaded Reference from User Code --------------- %
      NUMvod  =   ADIGATOR.NVAROFDIFF;
      fid     =   ADIGATOR.PRINT.FID;
      PFLAG   =   ADIGATOR.PRINT.FLAG;
      indent  =   ADIGATOR.PRINT.INDENT;
      if ADIGATOR.FORINFO.FLAG
        IncreaseForRefCount();
      end
      if ADIGATOR.EMPTYFLAG
        if length(s.subs) == 2
          y = cadaEmptyEval(x,s.subs{1},s.subs{2});
        else
          y = cadaEmptyEval(x,s.subs{1});
        end
        return
      elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
        y = ForSubsRef(x,s);
        return
      end
      xMrow = x.func.size(1); xNcol = x.func.size(2);
      % ------------------- Parse Reference Indices---------------------- %
      if length(s.subs) == 2
        [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1});
        [s.subs{2},subs2,logicflag2] = parseIndex(s.subs{2});
        if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
        if isempty(s.subs{2}) && isinf(xNcol); xNcol = 1; end
      else
        [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1});
        subs2 = []; logicflag2 = 0;
        if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
        if isempty(s.subs{1}) && isinf(xNcol); xNcol = 1; end
      end
      
      % ------------------Build Function Properties---------------------- %
      y.id = ADIGATOR.VARINFO.COUNT;
      
      if isinf(xMrow) 
        xvec=1; xnvec=2;
        ytemp = 1:xNcol;
        if length(s.subs)==2 && isequal(s.subs{1},':')
          ytemp = ytemp(s.subs{2});
          FMrow = Inf; FNcol = length(ytemp);
          snvec = 2; ynvec = 2;
        elseif length(s.subs)==1 && isequal(s.subs{1},':') && xNcol==1
          ytemp = 1; snvec =1;
          FMrow = Inf; FNcol=1;
          ynvec = 2;
        else
          error('Invalid vectorized subsref')
        end
      elseif isinf(xNcol); 
        xvec=2; xnvec=1;
        ytemp = (1:xMrow).';
        if length(s.subs)==2 && isequal(s.subs{2},':')
          ytemp = ytemp(s.subs{1});
          FMrow = length(ytemp); FNcol = Inf;
          ynvec = 1;
        elseif length(s.subs)==1 && isequal(s.subs{1},':') && xMrow==1
          ytemp = 1;
          FMrow = Inf; FNcol = 1;
          ynvec = 2;
        else
          error('Invalid vectorized subsref')
        end
        snvec = 1;
      else
        xvec=0;
        ytemp    = zeros(xMrow,xNcol);
        ytemp(:) = 1:xMrow*xNcol;
        ytemp    = ytemp(s.subs{:}); % Error here implies bad reference index
        [FMrow, FNcol] = size(ytemp);
      end
      isubs = ytemp(:);
      
      
      [funcstr,DPFLAG] = cadafuncname();
      y.func = struct('name',funcstr,'size',[FMrow,FNcol],'zerolocs',[],...
        'value',[]);
      if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 1
        AssgnForRefInds(isubs,s.subs);
      end
      
      % Function Numerics/Sparsity
      if ~isempty(x.func.value) && ~logicflag1 && ~logicflag2
        if xvec
          y.func.value = x.func.value(s.subs{snvec});
        else
          y.func.value = subsref(x.func.value,s);
        end
      elseif ~isempty(x.func.zerolocs) && ~logicflag1 && ~logicflag2
        if xvec
          xtemp = zeros(x.func.size(xnvec),1);
          xtemp(x.func.zerolocs) = 1;
          ytemp = xtemp(s.subs{xnvec});
          y.func.zerolocs = find(ytemp(:));
          if length(y.func.zerolocs) == y.func.size(ynvec)
            if ynvec ==1; y.func.value=zeros(FMrow,1); else y.func.value=zeros(1,FNcol); end
            y.func.zerolocs = [];
          end
        else
          xtemp = zeros(xMrow,xNcol);
          xtemp(x.func.zerolocs) = 1;
          ytemp = xtemp(s.subs{:});
          y.func.zerolocs = find(ytemp(:));
          if length(y.func.zerolocs) == FMrow*FNcol
            y.func.value = zeros(FMrow,FNcol);
            y.func.zerolocs = [];
          end
        end
      end
      if isfield(x.func,'logical')
        y.func.logical = [];
        if isempty(x.func.value) && (logicflag1 || logicflag2)
          error('Cannot use an unknown logical index to reference off of an unknown logical index');
        end
      end
      
      % Symbolic Logical Referencing
      if logicflag1 && logicflag2
        y.func.logicref = [subs1.id subs2.id];
      elseif logicflag1
        y.func.logicref = [subs1.id 0];
      elseif logicflag2
        y.func.logicref = [0 subs1.id];
      end
      
      
      % ------------------- Build Derivative Properties ----------------- %
      y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      for Vcount = 1:NUMvod;
        if ~isempty(x.deriv(Vcount).nzlocs)
          derivstr = cadadername(funcstr,Vcount);
          y.deriv(Vcount).name = derivstr;
          nzx = size(x.deriv(Vcount).nzlocs,1);
          if xvec
            dx = sparse(x.deriv(Vcount).nzlocs(:,1),x.deriv(Vcount).nzlocs(:,2),(1:nzx)',...
              x.func.size(xnvec),ADIGATOR.VAROFDIFF(Vcount).usize);
          else
            dx = sparse(x.deriv(Vcount).nzlocs(:,1),x.deriv(Vcount).nzlocs(:,2),(1:nzx)',...
              xMrow*xNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
          end

          [yrows,ycols,yind] = find(dx(isubs,:));
          if size(yrows,2) > 1
            yrows = yrows';
            ycols = ycols';
            yind = yind.';
          end
          y.deriv(Vcount).nzlocs = [yrows,ycols];
          if ~isempty(yind)
            if DPFLAG == 1
              % -------Derviatve Printing------
              Dind1 = cadaindprint(yind);
              if xvec
                fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(:,',Dind1,');\n']);
              else
                fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
              end
            end
          else
            y.deriv(Vcount).name = [];
          end
        end
      end
      
      % ---------------------- Function Printing ------------------------ %
      if PFLAG == 1
        if isempty(subs2)
          % Index referencing
          fprintf(fid,[indent,funcstr,' = ',x.func.name,'(',subs1.func.name,');\n']);
        else
          % Subs referencing
          fprintf(fid,[indent,funcstr,' = ',x.func.name,'(',subs1.func.name,',',subs2.func.name,');\n']);
        end
      end
      
%       % --------------- Check for Reference off of VOD ------------------ %
%       if isfield(x.func,'vodflag') && ...
%           isequal(y.deriv(x.func.vodflag).nzlocs(:,1),(1:FMrow*FNcol).')
%         y.func.vodflag = x.func.vodflag;
%       end
      % ---------------- Set Last Occurence and Var Count --------------- %
      if ~isempty(subs2)
        ADIGATOR.VARINFO.LASTOCC(subs2.id,1) = ADIGATOR.VARINFO.COUNT;
      end
      ADIGATOR.VARINFO.LASTOCC([subs1.id x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      y = class(y,'cada');
      if ADIGATOR.RUNFLAG == 1 && ADIGATOR.FORINFO.FLAG
        SubsrefUnion(x,y);
      end
    else % Reference from ADiGator module
      y =  y(s(scount).subs{:});
    end
  elseif strcmp(s(scount).type,'.') % Reference from ADiGator module
    if scount == 1
      if strcmp(s(scount).subs,'id')
        y = x.id;
      elseif strcmp(s(scount).subs,'func')
        y = x.func;
      elseif strcmp(s(scount).subs,'deriv')
        y = x.deriv;
      else
        error('invalid subscript reference for CADA')
      end
    else
      y = y.(s(scount).subs);
    end
  else
    error('invalid index reference for CADA')
  end
end
return
end

function [numeric, overloaded, logicflag] = parseIndex(index)
global ADIGATOR
if isa(index,'cada')
  overloaded = index;
  if ~isempty(index.func.value)
    numeric    = index.func.value;
    logicflag  = 0;
  elseif index.func.size(1) == 0 || index.func.size(2) == 0
    numeric    = [];
    logicflag  = 0;
  elseif isfield(index.func,'OnetoN')
    numeric = ':';
    logicflag  = 0;
  elseif isfield(index.func,'logical')
    overloaded = index;
    logicflag  = 1;
    if (isinf(index.func.size(1)) && index.func.size(2) == 1) || ...
        (isinf(index.func.size(2)) && index.func.size(1) == 1)
      if isequal(index.func.value,false)
        numeric = [];
        overloaded.func.name = '[]';
      else
        numeric = ':';
        overloaded.func.name = ':';
      end
    elseif any(isinf(index.func.size))
      error(['A vectorized logical reference index must be of dim ',...
        'N by 1 or 1 by N, where N is vectorized dim']);
    elseif ~isempty(index.func.zerolocs)
      numeric = true(index.func.size);
      numeric(index.func.zerolocs) = false;
      if ADIGATOR.PRINT.FLAG
        overloaded.func.name = cadaindprint(numeric);
      end
    else
      overloaded.func.name = ':';
      numeric = true(prod(overloaded.func.size),1);
    end
  else
    error('Cannot do strictly symbolic referencing/assignment.')
  end
elseif isnumeric(index)
  logicflag = 0;
  numeric = index;
  overloaded.id = [];
  if ADIGATOR.PRINT.FLAG
    overloaded.func.name = cadaindprint(index);
  end
elseif strcmp(index,':')
  logicflag = 0;
  numeric = index;
  overloaded.id = [];
  overloaded.func.name = ':';
else
  error('Invalid reference index.')
end

end

function IncreaseForRefCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF =...
  ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF + 1;
return
end

function AssgnForRefInds(inds,subs)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;

if isempty(inds)
  inds = -1;
end
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).INDICES = inds(:);
else
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).INDICES...
    (1:length(inds),ITERCOUNT) = inds;
end
% FLAGS - have one for subs1, one for subs2
%   1 => regular referencing
%   2 => logical referencing
%   3 => mixed between the two - NOT ALLOWED IF Y CHANGES SIZE
%   0 => if subs2 = 0, then it is linear index


if length(subs) == 2
  if isa(subs{1},'cada'); subs1 = subs{1}.func.value; else subs1 = subs{1}; end
  if islogical(subs1); Flag1 = 2; else Flag1 = 1; end
  if isa(subs{2},'cada'); subs2 = subs{2}.func.value; else subs2 = subs{2}; end
  if islogical(subs2); Flag2 = 2; else Flag2 = 1; end
else
  if isa(subs{1},'cada'); subs1 = subs{1}.func.name; else subs1 = subs{1}; end
  if islogical(subs1); Flag1 = 2; else Flag1 = 1; end
  Flag2 = 0;
end
if isempty(ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS)
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS = [0 0];
end
OFlag1 = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(1);
if OFlag1
  OFlag2 = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(2);
  if OFlag1 ~= Flag1
    ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(2) = 3;
  end
  if OFlag2 ~= Flag2
    ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(2) = 3;
  end
else
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS = [Flag1 Flag2];
end
return
end

function SubsrefUnion(x,y)
% Union y and assign x OPcount if hasnt been done yet.
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF;

if isempty(ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS)
  % Get overmapped x - we are probably already storing it somewhere, if we
  % arent, then store it as is because it wont be changing.
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS     = cell(1,2);
  if isempty(x.id)
    ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1} = [];
  else
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1} = x;
    end
  end
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{2}  = y;
else
  ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{2}  =...
    cadaUnionVars(ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{2},y);
  
  xOver = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1};
  if ~isempty(x.id)
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    if xOverLoc1 && xOver ~= xOverLoc1
      ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{1} = xOverLoc1;
    end
  end
end
ysize = y.func.size;
if any(ysize == 0) 
  ysize = [-1 -1];
end

Rsizes = [ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES,[ysize.';x.func.size.']];
ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES = Rsizes;
if (ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(1) == 3 || ...
    ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(2) == 3) && ...
    (any(Rsizes(1,1)~=Rsizes(1,:)) || any(Rsizes(2,1)~=Rsizes(2,:)))
  error(['References which 1. Change size on loop arrays, and 2. switch between',...
    'logical and index references, are not currently allowed.'])
end
return
end

function y = ForSubsRef(x,s)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF;
NUMvod    = ADIGATOR.NVAROFDIFF;
fid       = ADIGATOR.PRINT.FID;
PFLAG     = ADIGATOR.PRINT.FLAG;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0f',ADIGATOR.DERNUMBER);

y = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).VARS{2};
if isinf(y.func.size(1))
  yvec = 1;
  vecDim = ['size(',x.func.name,',1)'];
elseif isinf(y.func.size(2))
  vecDim = ['size(',x.func.name,',2)'];
  yvec = 2;
else
  yvec = 0;
end

[funcname,DPFLAG] = cadafuncname();
y.func.name = funcname;
y.id = ADIGATOR.VARINFO.COUNT;
CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
% ---------------------Derivative Printing------------------------------- %
if DPFLAG == 1
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      derivname = cadadername(funcname,Vcount);
      y.deriv(Vcount).name = derivname;
      if DPFLAG
        IndName  = ...
          ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).INDICES{Vcount,1};
        IndFlags = ...
          ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).INDICES{Vcount,3};
        nzy = size(y.deriv(Vcount).nzlocs,1);
        if yvec
          % Vectorized
          if IndFlags(1)
            % Indices change on this loop
            if IndFlags(2)
              % Indices are Sparse
              TD1 = ['cada',NDstr,'td1'];
              fprintf(fid,[indent,TD1,' = zeros(',vecDim,',%1.0d);\n'],nzy);
              fprintf(fid,[indent,TD1,'(:,logical(',IndName,'(:,',CountName,'))) = ',...
                x.deriv(Vcount).name,'(:,nonzeros(',IndName,'(:,',CountName,')));\n']);
              fprintf(fid,[indent,derivname,' = ',TD1,';\n']);
            else
              % Indices are non-Sparse
              fprintf(fid,[indent,derivname,' = ',...
                x.deriv(Vcount).name,'(',IndName,'(:,',CountName,'));\n']);
            end
          else
            % Indices do not change on this loop
            if IndFlags(2)
              % Indices are Sparse
              TD1 = ['cada',NDstr,'td1'];
              fprintf(fid,[indent,TD1,' = zeros(',vecDim,',1),%1.0d);\n'],nzy);
              fprintf(fid,[indent,TD1,'(:,logical(',IndName,')) = '...
                x.deriv(Vcount).name,'(:,nonzeros(',IndName,'));\n']);
              fprintf(fid,[indent,derivname, ' = ',TD1,';\n']);
            else
              % Indices are non-Sparse
              fprintf(fid,[indent,derivname,' = ',x.deriv(Vcount).name,'(:,',IndName,');\n']);
            end
          end
        else
          % Non-Vectorized
          if IndFlags(1)
            % Indices change on this loop
            if IndFlags(2)
              % Indices are Sparse
              TD1 = ['cada',NDstr,'td1'];
              fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],nzy);
              fprintf(fid,[indent,TD1,'(logical(',IndName,'(:,',CountName,'))) = ',...
                x.deriv(Vcount).name,'(nonzeros(',IndName,'(:,',CountName,')));\n']);
              fprintf(fid,[indent,derivname,' = ',TD1,';\n']);
            else
              % Indices are non-Sparse
              fprintf(fid,[indent,derivname,' = ',...
                x.deriv(Vcount).name,'(',IndName,'(:,',CountName,'));\n']);
            end
          else
            % Indices do not change on this loop
            if IndFlags(2)
              % Indices are Sparse
              TD1 = ['cada',NDstr,'td1'];
              fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],nzy);
              fprintf(fid,[indent,TD1,'(logical(',IndName,')) = '...
                x.deriv(Vcount).name,'(nonzeros(',IndName,'));\n']);
              fprintf(fid,[indent,derivname, ' = ',TD1,';\n']);
            else
              % Indices are non-Sparse
              fprintf(fid,[indent,derivname,' = ',x.deriv(Vcount).name,'(',IndName,');\n']);
            end
          end
        end
      end
    end
  end
end

% -----------------Parse Reference Inputs--------------------------- %
if length(s.subs) == 2
  if isa(s.subs{1},'cada')
    subs1 = s.subs{1};
  elseif isnumeric(s.subs{1})
    TFind1 = cadaindprint(s.subs{1});
    subs1.id = [];
    subs1.func = struct('name',TFind1,'size',size(s.subs{1}),...
      'zerolocs',[],'value',s.subs{1});
  elseif strcmp(s.subs{1},':')
    subs1.id = [];
    subs1.func = struct('name',':','size',[y.func.size(1), 1],...
      'zerolocs',[],'value',[]);
  else
    error('??? Invalid Reference')
  end
  if isa(s.subs{2},'cada')
    subs2 = s.subs{2};
  elseif isnumeric(s.subs{2})
    subs2.id = [];
    TFind1 = cadaindprint(s.subs{2});
    subs2.func = struct('name',TFind1,'size',size(s.subs{2}),...
      'zerolocs',[],'value',s.subs{2});
  elseif strcmp(s.subs{2},':')
    subs2.id = [];
    subs2.func = struct('name',':','size',[y.func.size(2), 1],...
      'zerolocs',[],'value',[]);
  else
    error('??? Invalid Reference')
  end
else
  if isa(s.subs{1},'cada')
    subs1 = s.subs{1};
  elseif isnumeric(s.subs{1})
    subs1.id = [];
    TFind1 = cadaindprint(s.subs{1});
    subs1.func = struct('name',TFind1,'size',size(s.subs{1}),...
      'zerolocs',[],'value',s.subs{1});
  elseif strcmp(s.subs{1},':')
    subs1.id = [];
    subs1.func = struct('name',':','size',[y.func.size(1), 1],...
      'zerolocs',[],'value',[]);
  else
    error('??? Invalid Reference')
  end
  subs2 = [];
end

% -------------------------Function Printing -----------------------------%
if PFLAG == 1
  if isempty(subs2) && ~isempty(ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{1,1})
    % Linear indexing off of a matrix which changes row size - need to fix
    % indexing - we dont care if columns change size.
    RowIndName = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{1,1}; % Rows change size
    RowDepFlag = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{1,3}(1); % Rows change size on this loop
    TFind1 = ['cada',NDstr,'tfind1']; TFind2 = ['cada',NDstr,'tfind2']; TFind3 = ['cada',NDstr,'tfind3'];
    fprintf(fid,[indent,TFind2,' = zeros(%1.0d,%1.0d);\n'],x.func.size(1),x.func.size(2));
    fprintf(fid,[indent,TFind2,'(:) = 1:%1.0d;\n'],x.func.size(1)*x.func.size(2));
    if RowDepFlag
      fprintf(fid,[indent,TFind3,' = ',TFind2,'(1:',RowIndName,'(',CountName,'),:);\n']);
    else
      fprintf(fid,[indent,TFind3,' = ',TFind2,'(1:',RowIndName,',:);\n']);
    end
    fprintf(fid,[indent,TFind1,' = ',TFind3,'(',subs1.func.name,');\n']);
    subs1.func.name = TFind1;
    % Set Logic Flag so that next section doesnt mess with the referencing
    % we just did.
    ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(1) = 2;
  end

  RowIndName = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{2,1};
  ColIndName = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{3,1};
  if ~isempty(RowIndName)
    RowDepFlag = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{2,3}(1);
  end
  if ~isempty(ColIndName)
    ColDepFlag = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).SIZES{3,3}(1);
  end
  if isempty(RowIndName) && isempty(ColIndName)
    % Neither Row or Column Size of Y Changes
    if isempty(subs2)
      % Index referencing
      fprintf(fid,[indent,y.func.name,' = ',x.func.name,'(',subs1.func.name,');\n']);
    else
      % Subs referencing
      fprintf(fid,[indent,y.func.name,' = ',x.func.name,'(',subs1.func.name,',',subs2.func.name,');\n']);
    end
  else
    RowLogicFlag = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(1);
    ColLogicFlag = ADIGATORFORDATA(INNERLOC).SUBSREF(REFCOUNT).FLAGS(2);
    if ~isempty(RowIndName) && RowDepFlag
      RowIndName = [RowIndName,'(',CountName,')'];
    end
    if ~isempty(ColIndName) && ColDepFlag
      ColIndName = [ColIndName,'(',CountName,')'];
    end
    tempfuncstr  = ['cada',NDstr,'tempf1'];

    if ~isempty(RowIndName) && ~isempty(ColIndName)
      % Rows and Columns Changing Size
      if isempty(subs2)
        % Single Index
        if RowLogicFlag == 2% Logical
          RHSstr = [x.func.name,'(',subs1.func.name,')'];
        elseif subs1.func.size(1) == 1 || subs1.func.size(2) == 1
          RHSstr = [x.func.name,'(',subs1.func.name,'(1:',RowIndName,'*',ColIndName,'))'];
        else
          RHSstr = [x.func.name,'(',subs1.func.name,'(1:',RowIndName,',1:',ColIndName,'))'];
        end
      else
        % Two Indices
        if strcmp(subs1.func.name,':')
          Ref1 = ['1:',RowIndName];
        elseif RowLogicFlag == 2
          Ref1 = subs1.func.name;
        else
          Ref1 = [subs1.func.name,'(1:',RowIndName,')'];
        end
        if strcmp(subs2.func.name,':')
          Ref2 = ['1:',ColIndName];
        elseif ColLogicFlag == 2
          Ref2 = subs2.func.name;
        else
          Ref2 = [subs2.func.name,'(1:',ColIndName,')'];
        end
        RHSstr = [x.func.name,'(',Ref1,',',Ref2,')'];
      end
      fprintf(fid,[indent,tempfuncstr,' = ',RHSstr,';\n']);
      if yvec == 1
        fprintf(fid,[indent,y.func.name,' = zeros(',vecDim,',%1.0f);\n'],y.func.size(2));
      elseif yvec == 2
        fprintf(fid,[indent,y.func.name,' = zeros(%1.0f,',vecDim,');\n'],y.func.size(1));
      else
        fprintf(fid,[indent,y.func.name,' = zeros(%1.0f,%1.0f);\n'],y.func.size(1),y.func.size(2));
      end
      fprintf(fid,[indent,y.func.name,'(1:',RowIndName,',1:',ColIndName,') = ',tempfuncstr,';\n']);
    elseif ~isempty(RowIndName)
      % Rows Changing Size
      if isempty(subs2)
        % Single Index
        if strcmp(subs1.func.name,':')
          RHSstr = [x.func.name,'(1:',RowIndName,')'];
        elseif RowLogicFlag == 2 || ColLogicFlag == 2
          RHSstr = [x.func.name,'(',subs1.func.name,')'];
        elseif subs1.func.size(2) > 1 && subs1.func.size(1) > 1
          RHSstr = [x.func.name,'(',subs1.func.name,'(1:',RowIndName,',:))'];
        else
          RHSstr = [x.func.name,'(',subs1.func.name,'(1:',RowIndName,'))'];
        end
      else
        % Two Indices
        if strcmp(subs1.func.name,':')
          Ref1 = ['1:',RowIndName];
        elseif RowLogicFlag == 2
          Ref1 = subs1.func.name;
        else
          Ref1 = [subs1.func.name,'(1:',RowIndName,')'];
        end
        RHSstr = [x.func.name,'(',Ref1,',',subs2.func.name,')'];
      end
      fprintf(fid,[indent,tempfuncstr,' = ',RHSstr,';\n']);
      if yvec
        fprintf(fid,[indent,y.func.name,' = zeros(%1.0f,',vecDim,');\n'],y.func.size(1));
      else
        fprintf(fid,[indent,y.func.name,' = zeros(%1.0f,%1.0f);\n'],y.func.size(1),y.func.size(2));
      end
      if y.func.size(2) > 1
        fprintf(fid,[indent,y.func.name,'(1:',RowIndName,',:) = ',tempfuncstr,';\n']);
      else
        fprintf(fid,[indent,y.func.name,'(1:',RowIndName,') = ',tempfuncstr,';\n']);
      end
    else
      % Columns Changing Size
      if isempty(subs2)
        % Single Index
        if RowLogicFlag == 2 || ColLogicFlag == 2
          RHSstr = [x.func.name,'(',subs1.func.name,')'];
        elseif subs1.func.size(2) > 1 && subs1.func.size(1) > 1
          RHSstr = [x.func.name,'(',subs1.func.name,'(:,1:',ColIndName,'))'];
        else
          RHSstr = [x.func.name,'(',subs1.func.name,'(1:',ColIndName,'))'];
        end
      else
        % Two Indices
        if strcmp(subs2.func.name,':')
          Ref2 = ['1:',ColIndName];
        elseif ColLogicFlag == 2
          Ref2 = subs2.func.name;
        else
          Ref2 = [subs2.func.name,'(1:',ColIndName,')'];
        end
        RHSstr = [x.func.name,'(',subs1.func.name,',',Ref2,')'];
      end
      fprintf(fid,[indent,tempfuncstr,' = ',RHSstr,';\n']);
      if yvec
        fprintf(fid,[indent,y.func.name,' = zeros(',vecDim,',%1.0f);\n'],y.func.size(2));
      else
        fprintf(fid,[indent,y.func.name,' = zeros(%1.0f,%1.0f);\n'],y.func.size(1),y.func.size(2));
      end
      if y.func.size(1) > 1
        fprintf(fid,[indent,y.func.name,'(:,1:',ColIndName,') = ',tempfuncstr,';\n']);
      else
        fprintf(fid,[indent,y.func.name,'(1:',ColIndName,') = ',tempfuncstr,';\n']);
      end
    end
  end
end

if ~isempty(subs2)
  ADIGATOR.VARINFO.LASTOCC(subs2.id,1) = ADIGATOR.VARINFO.COUNT;
end
ADIGATOR.VARINFO.LASTOCC([subs1.id x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
  y = class(y,'cada');
end

return
end