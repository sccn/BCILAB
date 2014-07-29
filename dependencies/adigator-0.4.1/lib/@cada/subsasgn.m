function y = subsasgn(x,s,b)
% CADA overloaded version of MATLAB function SUBASGN, to be used in
% conjuction with the ADiGator algorithm.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
ssize = length(s);
bscalarflag = 0;

if ssize == 1 && strcmp(s(1).type,'()')
  if ADIGATOR.FORINFO.FLAG
    IncreaseForAsgnCount();
  end
  if ADIGATOR.EMPTYFLAG
    if length(s.subs) == 2
      y = cadaEmptyEval(x,s.subs{1},s.subs{2},b);
    else
      y = cadaEmptyEval(x,s.subs{1},b);
    end
    return
  elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    y = ForSubsAsgn(x,s,b);
    return
  end
  % -----------------------Parse B/X input--------------------------------%
  if isnumeric(b)
    % b is numeric input
    [bMrow,bNcol] = size(b);
     btemp.id = [];
     btemp.func = struct('name',[],'size',[bMrow,bNcol],'zerolocs',[],...
       'value',b);
    if PFLAG
      if bMrow*bNcol == 1
        btemp.func.name = num2str(b,16);
      else
        cadamatprint(b,['cada',NDstr,'temp1']);
        btemp.func.name = ['cada',NDstr,'temp1'];
      end
    end
    btemp.func.size = [bMrow, bNcol];
    btemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    b = btemp;
    b = class(b,'cada');
  else
    bMrow = b.func.size(1);
    bNcol = b.func.size(2);
  end
  if isnumeric(x)
    % x is numeric input
    [xMrow,xNcol] = size(x);
     xtemp.id = [];
     xtemp.func = struct('name',[],'size',[xMrow,xNcol],'zerolocs',[],...
       'value',x);
    if PFLAG
      if xMrow*xNcol == 1
        xtemp.func.name = num2str(x,16);
      else
        cadamatprint(x,['cada',NDstr,'temp1']);
        xtemp.func.name = ['cada',NDstr,'temp1'];
      end
    end
    xtemp.func.size = [xMrow, xNcol];
    xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    x = xtemp;
    x = class(x,'cada');
  else
    xMrow = x.func.size(1);
    xNcol = x.func.size(2);
  end

  % ----------------------- Parse SUBS inputs --------------------------- %
  if length(s.subs) == 2
    [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1});
    [s.subs{2},subs2,logicflag2] = parseIndex(s.subs{2});
    if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
    if isempty(s.subs{2}) && isinf(xNcol); xNcol = 1; end
  else
    [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1});
    logicflag2 = 0; subs2 = [];
    if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
    if isempty(s.subs{1}) && isinf(xNcol); xNcol = 1; end
  end
  
  % ------------------ Check For Logical Assignment --------------------- %
  %logicflag1 => 1st assignment index is logical with unknown values
  %logicflag2 => 2nd assignment index is logical with unknown values
  if (logicflag1 || logicflag2) && (~bMrow || ~bNcol)
    error(['Cannot use logical assignment with ',...
      'unknown logical values to remove elements of an object'])
  end
  logicerrflag = 0;
  if bMrow*bNcol > 1
    % This could probably be done cleaner, but i think this covers all the
    % cases which are not allowed.
    if logicflag1
      if  ~isfield(b.func,'logicref')
        logicerrflag = 1;
      elseif b.func.logicref(1) ~= subs1.id
        if ~isempty(subs2) && prod(subs2.func.size) > 1
          logicerrflag = 1;
        elseif b.func.size(1) ~= 1 || b.func.logicref(2) ~= subs1.id
          logicerrflag = 1;
        end
      end
    end
    if logicflag2
      if  ~isfield(b.func,'logicref')
        logicerrflag = 1;
      elseif b.func.logicref(2) ~= subs2.id
        if prod(subs1.func.size) > 1
          logicerrflag = 1;
        elseif b.func.size(2) ~= 1 || b.func.logicref(1) ~= subs1.id
          logicerrflag = 1;
        end
      end
    end
  end
  if logicerrflag
    error(sprintf(['May only use logical assignment with unknown logical ',...
      'values if the object being assigned is either a scalar, or was ',...
      'created by a logical reference, where the logical reference ',...
      'index is the same as the logical assignment index.\n Example:\n',...
      'ind = x < a\nx(ind) = b(ind); is valid.\nx(x<a) = b(x<a); is not.'])); %#ok<SPERR>
  end
  %---------------------------------------------------------------------%
  %                      Build Function Properties                      %
  %---------------------------------------------------------------------%
  y.id = ADIGATOR.VARINFO.COUNT;

  
  if isinf(xMrow)
    ytemp = 1:xNcol;
    xvec = 1; xnvec = 2; ynvec = 2;
    vecDim = ['size(',x.func.name,',1)'];
    if isinf(bMrow) && isequal(s.subs{1},':')
      bnvec = 2; btemp = 1:bNcol;
    elseif isinf(bNcol) && bMrow==1 && isequal(s.subs{1},':')
      bnvec = 1; btemp = 1;
    elseif bMrow*bNcol==1 && ~cadaCheckForDerivs(b) && isequal(s.subs{1},':')
      bnvec = 0; btemp = 1;
    else
      error('Invalid vectorized subsasgn')
    end
    if length(s.subs) == 2
      ytemp(s.subs{2}) = btemp;
      snvec = 2;
    else
      snvec = 1;
    end
    FMrow = Inf; FNcol = length(ytemp);
  elseif isinf(xNcol)
    vecDim = ['size(',x.func.name,',2)'];
    ytemp = (1:xMrow).';
    xvec = 2; xnvec = 1; ynvec = 1;
    if isinf(bNcol) && ...
        (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
        (length(s.subs)==1 && isequal(s.subs{1},':'))
      bnvec = 1; btemp = (1:bMrow).';
    elseif isinf(bMrow) && bNcol==1 && ...
        (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
        (length(s.subs)==1 && isequal(s.subs{1},':'))
      bnvec = 2; btemp = 1;
    elseif bMrow*bNcol==1 && ~cadaCheckForDerivs(b) && ...
        (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
        (length(s.subs)==1 && isequal(s.subs{1},':'))
      bnvec = 0; btemp = 1;
    else
      error('Invalid vectorized subsasgn')
    end
    if length(s.subs) == 2
      ytemp(s.subs{1}) = btemp;
    end
    snvec = 1;
    FMrow = length(ytemp); FNcol = Inf;
  elseif prod(x.func.size) == 0 && isinf(bMrow)
    xvec = 1;
    vecDim = ['size(',b.func.name,',1)'];
    ytemp  = [];
    bnvec = 2; btemp = 1:bNcol;
    if length(s.subs) == 2 && isequal(s.subs{1},':')
      ynvec = 2; xnvec = 2;
      FMrow = Inf; FNcol = bNcol;
      snvec = 2;
    elseif length(s.subs) == 2 && isequal(s.subs{2},':') && bNcol == 1
      ynvec = 1; xnvec = 1;
      FMrow = 1; FNcol = Inf;
      snvec = 1;
    else
      error('Invalid vectorized subsasgn')
    end
    ytemp(s.subs{:}) = btemp;
  elseif prod(x.func.size) == 0 && isinf(bNcol)
    xvec = 1;
    vecDim = ['size(',b.func.name,',2)'];
    ytemp = [];
    bnvec = 1; btemp = 1:bMrow;
    if length(s.subs) == 2 && isequal(s.subs{1},':') && bMrow == 1
      xnvec = 2; ynvec = 2;
      FMrow = Inf; FNcol = 1;
      snvec = 2;
    elseif length(s.subs) == 2 && isequal(s.subs{2},':')
      xnvec = 1; ynvec = 1;
      FMrow = bMrow; FNcol = Inf;
      snvec = 1;
    else
      error('Invalid vectorized subsasgn')
    end
    ytemp(s.subs{:}) = btemp;
  else
    xvec = 0; bnvec =0;
    ytemp = zeros(x.func.size);
    btemp = zeros(b.func.size);
    ytemp(s.subs{:}) = btemp;
    [FMrow, FNcol] = size(ytemp);
  end

  % Check b scalar case to see if need to REPMAT
  if bMrow*bNcol ==1
    btemp = ytemp(s.subs{:});
    [bMrow,bNcol] = size(btemp);
    if bMrow*bNcol > 1
      bscalarflag = 1;
      b.func.size = [bMrow, bNcol];
    elseif xvec
      bscalarflag = 1;
    end
  end
  
  if xvec == 1
    ytemp = 1:FNcol;
    if length(s.subs)==2
      ytemp = ytemp(s.subs{2});
    end
  elseif xvec == 2
    ytemp = (1:FMrow).';
    if length(s.subs)==2
      ytemp = ytemp(s.subs{1});
    end
  else
    ytemp = zeros(FMrow,FNcol);
    ytemp(:) = (1:FMrow*FNcol)';
    ytemp = ytemp(s.subs{:});
  end
  isubs = ytemp(:);
  if xvec == 1
    insubs = 1:FNcol;
  elseif xvec==2
    insubs = 1:FMrow;
  else
    insubs = 1:FMrow*FNcol;
  end
  nsubs = length(isubs);
  insubs(isubs) = [];
  [funcstr,DPFLAG] = cadafuncname();
  y.func = struct('name',funcstr,'size',[FMrow,FNcol],'zerolocs',[],...
    'value',[]);
  if ADIGATOR.FORINFO.FLAG
    AssignForAsgnInds(isubs,s.subs,bscalarflag);
  end
  % -------------- Function Numeric Values and Sparsity ----------------- %
  if ~isempty(x.func.value) && ~isempty(b.func.value) &&...
      ~logicflag1 && ~logicflag2
    y.func.value = x.func.value;
    if xvec
      y.func.value(s.subs{snvec}) = b.func.value;
    else
      y.func.value(s.subs{:}) = b.func.value;
    end
  else
    spflag = 0;
    if ~isempty(x.func.value)
      xtemp = logical(x.func.value);
      spflag = 1;
    elseif ~isempty(x.func.zerolocs)
      if xvec
        xtemp = true(x.func.size(xnvec),1);
      else
        xtemp = true(xMrow,xNcol);
      end
      xtemp(x.func.zerolocs) = false;
      spflag = 1;
    elseif xvec
      xtemp = true(x.func.size(xnvec),1);
    else
      xtemp = true(xMrow,xNcol);
    end
    if ~isempty(b.func.value) && bscalarflag && ~b.func.value
      spflag = 1;
    elseif ~isempty(b.func.value) && ~bscalarflag
      btemp = logical(b.func.value);
      spflag = 1;
    elseif ~isempty(b.func.zerolocs) && ~bscalarflag
      if bnvec
        btemp = true(b.func.size(bnvec),1);
      else
        btemp = true(bMrow,bNcol);
      end
      btemp(b.func.zerolocs) = false;
      spflag = 1;
    elseif bnvec
      btemp = true(b.func.size(bnvec),1);
    else
      btemp = true(bMrow,bNcol);
    end
    
    if spflag || (~xvec && FMrow*FNcol > xMrow*xNcol) ||...
        (xvec && y.func.size(ynvec) > x.func.size(xnvec))
      ytemp = xtemp;
      if xvec
        ytemp(s.subs{snvec}) = btemp;
      else
        ytemp(s.subs{:}) = btemp;
      end
      if logicflag1 || logicflag2
        ytemp(isubs) = or(ytemp(isubs),reshape(btemp,size(ytemp(isubs))));
      end
      
      y.func.zerolocs = find(~ytemp(:));
      if ~xvec && length(y.func.zerolocs) == FMrow*FNcol
        y.func.zerolocs = [];
        y.func.value = zeros(FMrow,FNcol);
      elseif xvec && length(y.func.zerolocs) == y.func.size(ynvec)
        y.func.zerolocs = [];
        if xvec == 1
          y.func.value = zeros(1,FNcol);
        else
          y.func.value = zeros(FMrow,1);
        end
      end
    end
  end
  % NOTES FOR LOGICAL ASSIGNMENT: The way we returned the subs variables,
  % we placed a logical 1 where something MAY be true, but is not
  % necessarily true. Thus, the elements insubs are most definitely not
  % being changed, but the elements isubs MIGHT be changed, we dont know.
  if xvec == 1 && isinf(logicflag1)
    vecrefstr = subs1.func.name;
    if DPFLAG && cadaCheckForDerivs(x) && cadaCheckForDerivs(b)
      xvecrefstr = ['cada',NDstr,'tind2'];
      fprintf(fid,[indent,xvecrefstr,' = ~',vecrefstr,';\n']);
    end
  elseif xvec == 2 && isinf(logicflag2)
    vecrefstr = subs2.func.name;
    if DPFLAG && cadaCheckForDerivs(x) && cadaCheckForDerivs(b)
      xvecrefstr = ['cada',NDstr,'tind2'];
      fprintf(fid,[indent,xvecrefstr,' = ~',vecrefstr,';\n']);
    end
  elseif xvec
    vecrefstr = ':';
  end
  
  % ---------------------Build Derivative Properties--------------------- %
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  for Vcount = 1:NUMvod;
    if ~isempty(x.deriv(Vcount).nzlocs) && ~isempty(b.deriv(Vcount).nzlocs)
      % --------X AND B HAVE DERIVATIVES WRT VAR(VCOUNT)---------
      if bscalarflag && bMrow*bNcol > 1
        [b.deriv(Vcount).name,b.deriv(Vcount).nzlocs] = ...
          cadaRepDers(b.deriv(Vcount).name,b.deriv(Vcount).nzlocs,bMrow*bNcol,Vcount,DPFLAG);
      end
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      nzx = size(x.deriv(Vcount).nzlocs,1);
      nzb = size(b.deriv(Vcount).nzlocs,1);
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      xrows = x.deriv(Vcount).nzlocs(:,1);
      if xvec; ny = y.func.size(ynvec); else ny = FMrow*FNcol; end
      if ~xvec && FMrow > xMrow && xNcol > 1
        % X changes row size and is a matrix, need to re-order derivatives
        yref = zeros(FMrow,FNcol); yref(:) = 1:ny;
        xref = yref(1:xMrow,1:xNcol); xrows = xref(xrows);
      end
      xcols = x.deriv(Vcount).nzlocs(:,2);
      brows = b.deriv(Vcount).nzlocs(:,1);
      bcols = b.deriv(Vcount).nzlocs(:,2);
      brows = isubs(brows);
      if size(brows,2) > 1
        brows = brows';
      end
      dx = sparse(xrows,xcols,(1:nzx)',ny,nv);
      xindin = nonzeros(dx(insubs,:));
      if logicflag1 || logicflag2
        
        % Unknown logical assignment locations
        % dy how we have it now is actually fine - each element which is
        % being assigned may or may not be assigned, so dy how we have it
        % is the union of those two cases
        dy = sparse([xrows;brows],[xcols;bcols],ones(nzx+nzb,1),ny,nv);
        db = sparse(brows,bcols,(1:nzb)',ny,nv);
        nzy = nnz(dy);
        % dy is now the union of dx and db
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1
          yrows = yrows';
          ycols = ycols';
        end
        y.deriv(Vcount).nzlocs = [yrows,ycols];
        if DPFLAG
          % We have 3 possible cases here:
          % 1.) no vectorization - 1 or more asgn indices is logical
          % unknown
          %     - Project into matrix of size nsubs by nv
          % 2.) vectorization - non-vectorized index is logical unknown
          %     - Project into matrix of size N by nsubs*nv (N is vec dim)
          % 3.) vectorization - non-vectorized index is known, vectorized
          % index is logical unknown
          %     - Project into matrix of size N by nnz(dymat)
          if ~xvec
            mycase = 1;
          elseif logicflag1 == 1 || logicflag2 == 1
            mycase = 2;
          else
            mycase = 3;
          end
          dy  = sparse(yrows,ycols,1:nzy,ny,nv);
          dymat = dy(isubs,:);
          dxmat = dx(isubs,:);
          dbmat = db(isubs,:);
          [dxind1, ~, dxind2]  = find(dxmat(:));
          [dbind1, ~, dbind2]  = find(dbmat(:));
          [dyind1, ~, dyind2]  = find(dymat(:));
          % Set dymat to zeros
          dymatstr = ['cada',NDstr,'td1'];
          if mycase == 1
            fprintf(fid,[indent,dymatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
          elseif mycase == 2
            fprintf(fid,[indent,dymatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
          else
            fprintf(fid,[indent,dymatstr,' = zeros(',vecDim,',%1.0f);\n'],nnz(dymat));
            dymat2 = full(dymat);
            dymat2(logical(dymat)) = 1:nnz(dymat);
          end
          
          if any(any(dx(isubs,:)))
            % There are elements of dx which lie in the isubs row locs
            if mycase < 3
              lRefStr = ...
                getLogicalDerivRef(logicflag1,logicflag2,s,subs1,subs2,length(isubs),FMrow,FNcol,nv);
              dxasgnstr = ['cada',NDstr,'tind1'];
              fprintf(fid,[indent,dxasgnstr,' = ~',lRefStr,';\n']);
            end
            
            % Project dx into dxmat
            dxmatstr = ['cada',NDstr,'td2'];
            if mycase == 1
              fprintf(fid,[indent,dxmatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
            elseif mycase == 2
              fprintf(fid,[indent,dxmatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
            else
              fprintf(fid,[indent,dxmatstr,' = zeros(',vecDim,',%1.0f);\n'],nnz(dymat));
            end
            
            if ~isempty(xindin)
              % Not all elements of dx go into dxmat
              if mycase < 3
                dxind1str = cadaindprint(dxind1);
              else
                dxind1str = cadaindprint(full(dymat2(logical(dxmat))));
              end
              dxind2str = cadaindprint(dxind2);
              if mycase == 1
                fprintf(fid,[indent,dxmatstr,'(',dxind1str,') = ',x.deriv(Vcount).name,'(',dxind2str,');\n']);
              else
                fprintf(fid,[indent,dxmatstr,'(:,',dxind1str,') = ',x.deriv(Vcount).name,'(:,',dxind2str,');\n']);
              end
            else
              % All elements of dx go into dxmat
              dxind1 = dxind1(dxind2);
              if mycase < 3
                dxind1str = cadaindprint(dxind1);
              else
                dxind1str = cadaindprint(full(dymat2(logical(dxmat))));
              end
              if mycase == 1
                fprintf(fid,[indent,dxmatstr,'(',dxind1str,') = ',x.deriv(Vcount).name,';\n']);
              else
                fprintf(fid,[indent,dxmatstr,'(:,',dxind1str,') = ',x.deriv(Vcount).name,';\n']);
              end
            end
            
            % Assign elements of dxmat to dymat
            if mycase == 1
              fprintf(fid,[indent,dymatstr,'(',dxasgnstr,',:) = ',dxmatstr,'(',dxasgnstr,',:);\n']);
            elseif mycase == 2
              fprintf(fid,[indent,dymatstr,'(',xvecrefstr,',',dxasgnstr,') = ',dxmatstr,'(',xvecrefstr,',',dxasgnstr,');\n']);
            else
              fprintf(fid,[indent,dymatstr,'(',xvecrefstr,',:) = ',dxmatstr,'(',xvecrefstr,',:);\n']);
            end
          elseif mycase < 3
            lRefStr = ...
              getLogicalDerivRef(logicflag1,logicflag2,s,subs1,subs2,length(isubs),FMrow,FNcol,nv);
          end
          
          % Project db into dbmat and assign elements of dbmat to dymat
          dbmatstr = ['cada',NDstr,'td3'];
          
          dbind1 = dbind1(dbind2);
          if mycase < 3
            dbind1str = cadaindprint(dbind1);
            dbasgnstr = lRefStr;
          end
          if mycase == 1
            fprintf(fid,[indent,dbmatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
            fprintf(fid,[indent,dbmatstr,'(',dbind1str,') = ',b.deriv(Vcount).name,';\n']);
            fprintf(fid,[indent,dymatstr,'(',dbasgnstr,',:) = ',dbmatstr,'(',dbasgnstr,',:);\n']);
          elseif mycase == 2
            fprintf(fid,[indent,dbmatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
            fprintf(fid,[indent,dbmatstr,'(:,',dbind1str,') = ',b.deriv(Vcount).name,';\n']);
            fprintf(fid,[indent,dymatstr,'(',vecrefstr,',',dbasgnstr,') = ',dbmatstr,'(',vecrefstr,',',dbasgnstr,');\n']);
          else
            fprintf(fid,[indent,dbmatstr,' = zeros(',vecDim,',%1.0f);\n'],nnz(dbmat));
            dbind1str = cadaindprint(full(dymat2(logical(dbmat))));
            dbind2str = cadaindprint(dbind2);
            fprintf(fid,[indent,dbmatstr,'(:,',dbind1str,') = ',b.deriv(Vcount).name,'(:,',dbind2str,');\n']);
            fprintf(fid,[indent,dymatstr,'(',vecrefstr,',:) = ',dbmatstr,'(',vecrefstr,',:);\n']);
          end
          
          % Assign dy
          
          if ~isempty(xindin)
            % Assign elements of dx which are for sure not changed.
            DYstr = ['cada',NDstr,'td4'];
            if mycase == 1
              fprintf(fid,[indent,DYstr,' = zeros(%1.0f,1);\n'],nzy);
            else
              fprintf(fid,[indent,DYstr,' = zeros(',vecDim,',%1.0f);\n'],nzy);
            end
            xindout = nonzeros(dy(insubs,:));
            Dind1 = cadaindprint(xindout);
            Dind2 = cadaindprint(xindin);
            if mycase == 1
              fprintf(fid,[indent,DYstr,'(',Dind1,') = ',x.deriv(Vcount).name,'(',Dind2,');\n']);
            else
              fprintf(fid,[indent,DYstr,'(:,',Dind1,') = ',x.deriv(Vcount).name,'(:,',Dind2,');\n']);
            end
            
            % Now assign elements from dymatstr to dy
            dyind2str = cadaindprint(dyind2);
            if mycase < 3
              dyind1str = cadaindprint(dyind1);
              if mycase == 1
                fprintf(fid,[indent,DYstr,'(',dyind2str,') = ',dymatstr,'(',dyind1str,');\n']);
              else
                fprintf(fid,[indent,DYstr,'(:,',dyind2str,') = ',dymatstr,'(:,',dyind1str,');\n']);
              end
            else
              fprintf(fid,[indent,DYstr,'(:,',dyind2str,') = ',dymatstr,';\n']);
            end
            fprintf(fid,[indent,derivstr,' = ',DYstr,';\n']);
          elseif mycase < 3
            % Just get dy from dymat
            dyind1str = cadaindprint(dyind1);
            if mycase == 1
              fprintf(fid,[indent,derivstr,' = ',dymatstr,'(',dyind1str,');\n']);
            else
              fprintf(fid,[indent,derivstr,' = ',dymatstr,'(:,',dyind1str,');\n']);
            end
          else
            % dy is equal to dymat
            fprintf(fid,[indent,derivstr,' = ',dymatstr,';\n']);
          end
        end
      elseif ~isempty(xindin)
        dy = sparse([xrows;brows],[xcols;bcols],ones(nzx+nzb,1),ny,nv);
        % dy is now the union of dx and db
        % Regular assignment - know assignment locations
        diffref = nonzeros(dx(isubs,:));
        xdiffr = xrows(diffref);
        xdiffc = xcols(diffref);
        xdiff = sparse(xdiffr,xdiffc,ones(numel(xdiffr),1),ny,nv);
        % xdiff is
        dy = dy-xdiff;
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1
          yrows = yrows';
          ycols = ycols';
        end
        y.deriv(Vcount).nzlocs = [yrows,ycols];
        % ---------Derivative Printing-----------
        if DPFLAG == 1
          nzy = nnz(dy);
          dy = sparse(yrows,ycols,(1:nzy)',ny,nv);
          xindout = nonzeros(dy(insubs,:));
          %bindout = nonzeros(dy(isubs,:));
          db = sparse(brows,bcols,1:nzb,ny,nv);
          bindout = full(dy(logical(db)));
          bindin = nonzeros(db);
          bindinflag = ~isequal(bindin,(1:nzb).');
          TD1 = ['cada',NDstr,'td1'];
          Dind1 = cadaindprint(bindout);
          Dind2 = cadaindprint(xindout);
          Dind3 = cadaindprint(xindin);
          if bindinflag
            Dind4 = cadaindprint(bindin);
          end
          if xvec
            fprintf(fid,[indent,TD1,' = zeros(',vecDim,',%1.0f);\n'],nzy);
            if bindinflag
              fprintf(fid,[indent,TD1,'(:,',Dind1,') = ',b.deriv(Vcount).name,'(:,',Dind4,');\n']);
            else
              fprintf(fid,[indent,TD1,'(:,',Dind1,') = ',b.deriv(Vcount).name,';\n']);
            end
            fprintf(fid,[indent,TD1,'(:,',Dind2,') = ',x.deriv(Vcount).name,'(:,',Dind3,');\n']);
          else
            fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nzy);
            if bindinflag
              fprintf(fid,[indent,TD1,'(',Dind1,') = ',b.deriv(Vcount).name,'(',Dind4,');\n']);
            else
              fprintf(fid,[indent,TD1,'(',Dind1,') = ',b.deriv(Vcount).name,';\n']);
            end
            fprintf(fid,[indent,TD1,'(',Dind2,') = ',x.deriv(Vcount).name,'(',Dind3,');\n']);
          end
          fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
        end
      else
        dy = sparse(brows,bcols,ones(nzb,1),ny,nv);
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1
          yrows = yrows';
          ycols = ycols';
        end
        y.deriv(Vcount).nzlocs = [yrows,ycols];
        if DPFLAG == 1
          db = sparse(brows,bcols,1:nzb,ny,nv);
          bindin = nonzeros(db);
          if isequal(bindin,(1:nzb).');
            fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,';\n']);
          else
            Dind1 = cadaindprint(bindin);
            if xvec
              fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,'(:,',Dind1,');\n']);
            else
              fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,'(',Dind1,');\n']);
            end
          end
        end
      end
    elseif ~isempty(x.deriv(Vcount).nzlocs)
      % ---------------X has derivatives------------- %
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      nzx = size(x.deriv(Vcount).nzlocs,1);
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      if xvec; ny = y.func.size(ynvec); else ny = FMrow*FNcol; end
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2);
      if ~xvec && FMrow < xMrow && xNcol > 1
        % Removing rows from a matrix
        yref = zeros(FMrow,FNcol); yref(:) = 1:ny;
        yref(xMrow,xNcol) = 0; xnindex = 1:nzx;
        xrows = yref(xrows); xnindex = xnindex(logical(xrows));
        xcols = xcols(logical(xrows)); xrows = nonzeros(xrows);
        dx = sparse(xrows,xcols,xnindex,ny,nv);
        [xrows,xcols,xnindex] = find(dx);
        if size(xrows,2) > 1; xrows = xrows.'; xcols = xcols.'; end
        if ~isequal(xrows,x.deriv(Vcount).nzlocs(:,1))
          y.deriv(Vcount).nzlocs = [xrows,xcols];
          if DPflag
            Dind1 = cadaindprint(xnindex);
            fprintf(fid,[indent,y.deriv(Vcount).name,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
          end
        else
          y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
        end
      else
        % Assigning something with no derivatives.
        if ~xvec && FMrow > xMrow && xNcol > 1
          % X changes row size and is a matrix, need to re-order derivatives
          yref = zeros(FMrow,FNcol); yref(:) = 1:ny;
          xref = yref(1:xMrow,1:xNcol); xrows = xref(xrows);
          if size(xrows,2) > 1; xrows = xrows.'; end
        end
        dx = sparse(xrows,xcols,(1:nzx)',ny,nv);
        dxmat = dx(isubs,:);
        if nnz(dxmat)
          % Some derivatives got removed.
          if logicflag1 || logicflag2
            % Which derivatives get removed is determined by an unknown
            % logical index
            % In this case, dx as we have it is the sparsity pattern of dy,
            % though its possible that in some cases it could have
            % derivatives removed.
            % We have 3 possible cases here:
            % 1.) no vectorization - 1 or more asgn indices is logical
            % unknown
            %     - Project into matrix of size nsubs by nv
            % 2.) vectorization - non-vectorized index is logical unknown
            %     - Project into matrix of size N by nsubs*nv (N is vec dim)
            % 3.) vectorization - non-vectorized index is known, vectorized
            % index is logical unknown
            %     - Project into matrix of size N by nnz(dymat)
            if ~xvec
              mycase = 1;
            elseif logicflag1 == 1 || logicflag2 == 1
              mycase = 2;
            else
              mycase = 3;
            end
            [yrows,ycols] = find(dx);
            if size(yrows,2) > 1
              yrows = yrows';
              ycols = ycols';
            end
            y.deriv(Vcount).nzlocs = [yrows ycols];
            if DPFLAG
              if mycase < 3
                lRefStr = ...
                  getLogicalDerivRef(logicflag1,logicflag2,s,subs1,subs2,length(isubs),FMrow,FNcol,nv);
              end
              dy = sparse(yrows,ycols,(1:nzx),ny,nv);
              dymat    = dy(isubs,:);
              [dyind1, ~, dyind2]  = find(dymat(:));
              % Project elements of dx into dymat
              dymatstr = ['cada',NDstr,'td1'];
              if mycase == 1
                fprintf(fid,[indent,dymatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
              elseif mycase == 2
                fprintf(fid,[indent,dymatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
              end

              [dxind1, ~, dxind2]  = find(dxmat(:));
              if nnz(dx(insubs,:))
                % Some elements of dx lie in subs locs, some dont
                % Project dx into dymat and remove elements
                dxind2str = cadaindprint(dxind2);
                if mycase < 3
                  dxind1str = cadaindprint(dxind1);
                  if mycase == 1
                    fprintf(fid,[indent,dymatstr,'(',dxind1str,') = ',x.deriv(Vcount).name,'(',dxind2str,');\n']);
                    fprintf(fid,[indent,dymatstr,'(',lRefStr,',:) = 0;\n']);
                  else
                    fprintf(fid,[indent,dymatstr,'(:,',dxind1str,') = ',x.deriv(Vcount).name,'(:,',dxind2str,');\n']);
                    fprintf(fid,[indent,dymatstr,'(',vecrefstr,',',lRefStr,') = 0;\n']);
                  end
                else
                  fprintf(fid,[indent,dymatstr,' = ',x.deriv(Vcount).name,'(:,',dxind2str,');\n']);
                  fprintf(fid,[indent,dymatstr,'(',vecrefstr,',:) = 0;\n']);
                end

                % First assign elements of dx which are for sure not changed.
                DYstr = ['cada',NDstr,'td2'];
                if mycase == 1
                  fprintf(fid,[indent,DYstr,' = zeros(%1.0f,1);\n'],nzy);
                else
                  fprintf(fid,[indent,DYstr,' = zeros(',vecDim,',%1.0f);\n'],nzy);
                end
                xindout = nonzeros(dy(insubs,:));
                Dind1 = cadaindprint(xindout);
                Dind2 = cadaindprint(xindin);
                if mycase == 1
                  fprintf(fid,[indent,DYstr,'(',Dind1,') = ',x.deriv(Vcount).name,'(',Dind2,');\n']);
                else
                  fprintf(fid,[indent,DYstr,'(:,',Dind1,') = ',x.deriv(Vcount).name,'(:,',Dind2,');\n']);
                end

                % Extract derivs from dysubs
                dyind2str = cadaindprint(dyind2);
                if mycase == 1
                  dyind1str = cadaindprint(dyind1);
                  fprintf(fid,[indent,DYstr,'(',dyind2str,') = ',dymatstr,'(',dyind1str,');\n']);
                elseif mycase == 2
                  dyind1str = cadaindprint(dyind1);
                  fprintf(fid,[indent,DYstr,'(:,',dyind2str,') = ',dymatstr,'(:,',dyind1str,');\n']);
                else
                  fprintf(fid,[indent,DYstr,'(',dyind2str,') = ',dymatstr,';\n']);
                end
                fprintf(fid,[indent,derivstr,' = ',DYstr,';\n']);
              else
                % All elements of dx lie in subs locations
                if mycase < 3
                  dxind1 = dxind1(dxind2);
                  dxind1str = cadaindprint(dxind1);
                  dyind1str = cadaindprint(dyind1);
                  if mycase == 1
                    fprintf(fid,[indent,dymatstr,'(',dxind1str,') = ',x.deriv(Vcount).name,';\n']);
                    fprintf(fid,[indent,dymatstr,'(',lRefStr,',:) = 0;\n']);
                    fprintf(fid,[indent,derivstr,' = ',dymatstr,'(',dyind1str,');\n']);
                  else
                    fprintf(fid,[indent,dymatstr,'(:',dxind1str,') = ',x.deriv(Vcount).name,';\n']);
                    fprintf(fid,[indent,dymatstr,'(',vecrefstr,',',lRefStr,') = 0;\n']);
                    fprintf(fid,[indent,derivstr,' = ',dymatstr,'(',dyind1str,');\n']);
                  end
                else
                  dxind2str = cadaindprint(dxind2);
                  fprintf(fid,[indent,dymatstr,' = ',x.deriv(Vcount).name,'(:,',dxind2str,');\n']);
                  fprintf(fid,[indent,dymatstr,'(',vecrefstr,',:) = 0;\n']);
                  fprintf(fid,[indent,derivstr,' = ',dymatstr,';\n']);
                end
              end
            end
          else
            diffref = nonzeros(dx(isubs,:));
            xdiffr = xrows(diffref);
            xdiffc = xcols(diffref);
            dy = sparse([xrows;xdiffr],[xcols;xdiffc],[ones(nzx,1);-ones(numel(xdiffr),1)],ny,nv);
            xind = nonzeros(dx(insubs,:));
            if ~isempty(xind)
              [yrows,ycols] = find(dy);
              if size(yrows,2) > 1
                yrows = yrows';
                ycols = ycols';
              end
              y.deriv(Vcount).nzlocs = [yrows,ycols];
              if DPFLAG == 1
                Dind1 = cadaindprint(xind);
                if xvec
                  fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(:,',Dind1,');\n']);
                else
                  fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
                end
              end
            else
              % All derivatives got removed.
              y.deriv(Vcount).name = [];
            end
          end
        elseif ~isequal(x.deriv(Vcount).nzlocs(:,1),xrows)
            % Derivatives got re-ordered from the change in row size.
            [xrows,xcols,xnindex] = find(dx);
            if size(xrows,2) > 1; xrows = xrows.'; xcols = xcols.'; end
            y.deriv(Vcount).nzlocs = [xrows,xcols];
            if DPFLAG
              Dind1 = cadaindprint(xnindex);
              fprintf(fid,[indent,y.deriv(Vcount).name,' = ',x.deriv(Vcount).name,'(',Dind1,');\n']);
            end
        else
          y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
          % Dont need to print since y has same name as x.
        end
      end
    elseif ~isempty(b.deriv(Vcount).nzlocs)
      % B has derivatives
      if bscalarflag && bMrow*bNcol > 1
        [b.deriv(Vcount).name,b.deriv(Vcount).nzlocs] =...
          cadaRepDers(b.deriv(Vcount).name,b.deriv(Vcount).nzlocs,bMrow*bNcol,Vcount,DPFLAG);
      end
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      brows = b.deriv(Vcount).nzlocs(:,1);
      bcols = b.deriv(Vcount).nzlocs(:,2);
      nzb = length(brows);
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      if xvec; ny = y.func.size(ynvec); else ny = FMrow*FNcol; end
      brows = isubs(brows);
      if size(brows,2) > 1
        brows = brows';
      end
      db = sparse(brows,bcols,1:nzb,ny,nv);
      [yrows, ycols, bind] = find(db);
      if size(yrows,2) > 1
        yrows = yrows';
        ycols = ycols';
        bind  = bind.';
      end
      y.deriv(Vcount).nzlocs = [yrows,ycols];
      if DPFLAG == 1
        if logicflag1 || logicflag2
          % We have 3 possible cases here:
          % 1.) no vectorization - 1 or more asgn indices is logical
          % unknown
          %     - Project into matrix of size nsubs by nv
          % 2.) vectorization - non-vectorized index is logical unknown
          %     - Project into matrix of size N by nsubs*nv (N is vec dim)
          % 3.) vectorization - non-vectorized index is known, vectorized
          % index is logical unknown
          %     - Project into matrix of size N by nnz(dymat)
          if ~xvec
            mycase = 1;
          elseif logicflag1 == 1 || logicflag2 == 1
            mycase = 2;
          else
            mycase = 3;
          end
          if mycase < 3
            lRefStr = ...
              getLogicalDerivRef(logicflag1,logicflag2,s,subs1,subs2,length(isubs),FMrow,FNcol,nv);
          end
          dy = sparse(yrows,ycols,(1:nzb),ny,nv);
          dymat    = dy(isubs,:);
          dbmat    = db(isubs,:);
          [dyind1, ~, dyind2]  = find(dymat(:));
          % Make dymat and dbmat
          dymatstr = ['cada',NDstr,'td1'];
          dbmatstr = ['cada',NDstr,'td2'];
          if mycase == 1
            fprintf(fid,[indent,dymatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
            fprintf(fid,[indent,dbmatstr,' = zeros(%1.0f,%1.0f);\n'],nsubs,nv);
          elseif mycase == 2
            fprintf(fid,[indent,dymatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
            fprintf(fid,[indent,dbmatstr,' = zeros(',vecDim,',%1.0f);\n'],nsubs*nv);
          else
            fprintf(fid,[indent,dymatstr,' = zeros(',vecDim,',%1.0f);\n'],nnz(dymat));
          end
          
          % Project db into dbmat
          [dbind1, ~, dbind2]  = find(dbmat(:));
          
          if mycase < 3
            dbind1 = dbind1(dbind2);
            dbind1str = cadaindprint(dbind1);
            dbasgnstr = lRefStr;
            dyind1str = cadaindprint(dyind1);
            if mycase == 1
              fprintf(fid,[indent,dbmatstr,'(',dbind1str,') = ',b.deriv(Vcount).name,';\n']);
              fprintf(fid,[indent,dymatstr,'(',dbasgnstr,',:) = ',dbmatstr,'(',dbasgnstr,',:);\n']);
              fprintf(fid,[indent,derivstr,' = ',dymatstr,'(',dyind1str,');\n']);
            else
              fprintf(fid,[indent,dbmatstr,'(:,',dbind1str,') = ',b.deriv(Vcount).name,';\n']);
              fprintf(fid,[indent,dymatstr,'(',vecrefstr,',',dbasgnstr,') = ',dbmatstr,'(',vecrefstr,',',dbasgnstr,');\n']);
              fprintf(fid,[indent,derivstr,' = ',dymatstr,'(:,',dyind1str,');\n']);
            end
          else
            dbind2str = cadaindprint(dbind2);
            fprintf(fid,[indent,dbmatstr,' = ',b.deriv(Vcount).name,'(:,',dbind2str,');\n']);
            fprintf(fid,[indent,dymatstr,'(',vecrefstr,',:) = ',dbmatstr,'(',vecrefstr,',:);\n']);
            fprintf(fid,[indent,derivstr,' = ',dymatstr,';\n']);
          end
        elseif isequal(bind,(1:nzb).')
          fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,';\n']);
        else
         % Derivatives of b get reordered
          Dind1 = cadaindprint(bind);
          if xvec
            fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,'(:,',Dind1,');\n']);
          else
            fprintf(fid,[indent,derivstr,' = ',b.deriv(Vcount).name,'(',Dind1,');\n']);
          end
        end
      end
    end
  end
  
  %-------------------Function Printing---------------%
  if PFLAG == 1
    if ~isempty(x.func.name) && ~strcmp(x.func.name,y.func.name)
      fprintf(fid,[indent,funcstr,' = ',x.func.name,';\n']);
    end
    if logicflag1 && logicflag2 && ~bscalarflag
      fprintf(fid,[indent,funcstr,'(',subs1.func.name,',',subs2.func.name,')',...
        ' = ',b.func.name,'(',subs1.func.name,',',subs2.func.name,');\n']);
    elseif logicflag1 && ~bscalarflag
      if isempty(subs2)
        fprintf(fid,[indent,funcstr,'(',subs1.func.name,')',...
          ' = ',b.func.name,'(',subs1.func.name,');\n']);
      else
        fprintf(fid,[indent,funcstr,'(',subs1.func.name,',',subs2.func.name,')',...
          ' = ',b.func.name,'(',subs1.func.name,',:);\n']);
      end
    elseif logicflag2 && ~bscalarflag
      fprintf(fid,[indent,funcstr,'(',subs1.func.name,',',subs2.func.name,')',...
          ' = ',b.func.name,'(:,',subs1.func.name,');\n']);
    elseif size(s.subs,2) == 1
      fprintf(fid,[indent,funcstr,'(',subs1.func.name,') = ',b.func.name,';\n']);
    else
      fprintf(fid,[indent,funcstr,'(',subs1.func.name,',',subs2.func.name,') = ',b.func.name,';\n']);
    end
  end
  if ~isempty(subs2)
    ADIGATOR.VARINFO.LASTOCC(subs2.id,1) = ADIGATOR.VARINFO.COUNT;
  end
  ADIGATOR.VARINFO.LASTOCC([subs1.id y.id x.id b.id],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');

  if ADIGATOR.FORINFO.FLAG
    SubsasgnUnion(y,b);
  end
elseif strcmp(s(1).type,'.')
  if strcmp(s(1).subs,'id')
    y = x; y.id = b;
  elseif strcmp(s(1).subs,'func')
    y = x;
    if ssize == 1
      y.func = b;
    else
      y.func = subsasgn(x.func,s(2:end),b);
    end
  elseif strcmp(s(1).subs,'deriv')
    y = x;
    if ssize == 1
      y.deriv = b;
    else
      y.deriv = subsasgn(x.deriv,s(2:end),b);
    end
  elseif ~prod(x.func.size)
    y = [];
    y.(s(1).subs) = b;
  else
    error('Inncorrect subsasgn method for type cada');
  end
else
  error('Inncorrect subsasgn method for type cada');
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
        logicflag = 0;
      else
        numeric = ':';
        logicflag = Inf;
      end
    elseif any(isinf(index.func.size))
      error(['A vectorized logical reference index must be of dim ',...
        'N by 1 or 1 by N, where N is vectorized dim']);
    elseif ~isempty(index.func.zerolocs)
      numeric = true(index.func.size);
      numeric(index.func.zerolocs) = false;
    else
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

function IncreaseForAsgnCount()
% Increase the ASGN count
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN =...
  ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN + 1;
return
end

function AssignForAsgnInds(inds,subs,scalarflag)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN;
if isempty(inds) 
  inds = inf;
end
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES = inds(:);
else
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES...
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
  if isa(subs{1},'cada'); subs1 = subs{1}.func.value; else subs1 = subs{1}; end
  if islogical(subs1); Flag1 = 2; else Flag1 = 1; end
  Flag2 = 0;
end

if isempty(ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS)
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS = [0 0 0];
end
OFlag1 = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(1);
if length(inds) == 1
  scalarflag = 2;
end
if OFlag1
  OFlag2 = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(2);
  if OFlag1 ~= Flag1
    ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(2) = 3;
  end
  if OFlag2 ~= Flag2
    ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(2) = 3;
  end
  
  
  if scalarflag ~= ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(3)
    if scalarflag == 2
      % Its cool, leave the flag alone.
    elseif ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(3) == 2
      ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(3) = scalarflag;
    else
      error(sprintf(['Currently cannot do SUBSASGN in a for loop, when the variable\n',...
        'coming in changes from a scalar to a vector, with the same assignment index.\n',...
        'Ex. for I = 1:N\n    if I < 3; b = x(1); else b = x(1:N); end\n    y(1:N) = b\n',...
        'This statement would not be allowed, but instead could write\n',...
        'if I== 1; b = repmat(x(1),N,1); else b = x(1:N);'])) %#ok<SPERR>
    end
  end
else
  % First assignment
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS = [Flag1 Flag2 scalarflag];
end
return
end

function SubsasgnUnion(x,b)
% create b union, assign x if it hasn't been already.
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN;
if ~isa(b,'cada'); b = class(b,'cada'); end
if isempty(ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS)
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS = cell(1,2);
%  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{1} = ...
%    ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
  % Can see if we are storing b already or not.
  if ~isempty(b.id) && ADIGATOR.VARINFO.NAMELOCS(b.id,1)
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    bOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(b.id,1);
    bOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(b.id,2);
    if bOverLoc1 && b.id >= StartCount && b.id <= EndCount
      ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} = bOverLoc1;
    elseif bOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==bOverLoc2)
      ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} = bOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} = b;
    end
  else
    ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} = b;
  end
elseif isa(ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2},'cada')
  % Is an intermediate variable coming into this, we have no idea how it
  % will change, so just change it here.
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} =...
    cadaUnionVars(b,ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2});
else
  bOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(b.id,1);
  bOver     = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2};
  if bOverLoc1 && bOver ~= bOverLoc1
    ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2} = bOverLoc1;
  end
end
ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{1} = ...
    ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
x.func.size(isinf(x.func.size)) = 1;
b.func.size(isinf(b.func.size)) = 1;

bsize = b.func.size;
bsize(~logical(bsize)) = inf;
if ~isempty(ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES)
  bsizes = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES(3,:);
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES = ...
    [ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES,[x.func.size.';bsize.']];
  if (ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(1) == 3 || ...
      ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(2) == 3) && ...
      (any(bsizes(3,1) ~= bsizes(3,:)) || ...
      any(bsizes(3,1) ~= bsizes(3,:)))
    error(['References which 1. Change size on loop arrays, and 2. switch between',...
      'logical and index references, are not currently allowed.'])
  end
else
  ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES = [x.func.size.';bsize.'];
end
return
end

function y = ForSubsAsgn(x,s,b)
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN;
NUMvod    = ADIGATOR.NVAROFDIFF;
fid       = ADIGATOR.PRINT.FID;
PFLAG     = ADIGATOR.PRINT.FLAG;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% -------------------Get y from GLOBALFORDATA---------------------------- %
xOverLoc = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{1};
xOver    = ADIGATORVARIABLESTORAGE.OVERMAP{xOverLoc};
x = cadaPrintReMap(x,xOver,x.id);

y.id    = ADIGATOR.VARINFO.COUNT;
y.func  = x.func;
y.deriv = x.deriv;
[funcstr,DPFLAG]  = cadafuncname();
y.func.name       = funcstr;

% -----------------Parse SUBS inputs---------------------------
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
% -----------------------Parse B input--------------------------------%
if isnumeric(b)
  % b is numeric input
  [bMrow,bNcol] = size(b);
  btemp.id = [];
  if bMrow*bNcol == 1
    btemp.func.name = num2str(b,16);
  else
    btemp.func.name = cadamatprint(b);
  end
  btemp.func.size   = [bMrow, bNcol];
  btemp.deriv       = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  btemp.func.value  = b;
  btemp.func.zerolocs = [];
  b = btemp;
  b = class(b,'cada');
else
  bMrow = b.func.size(1);
  bNcol = b.func.size(2);
end

% ---Get OverMapped B--- %
bOver = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).VARS{2};

bscalarflag = 0;
% Check for the bscalarflag case.
if bMrow*bNcol == 1 && bOver.func.size(1)*bOver.func.size(2) > 1
  % The input b is scalar, and the overmapped is not. Either the OverMap
  % went bad, or this is bscalarflag case.
  if isempty(subs2)
    if subs1.func.size(1) > 1 || subs1.func.size(2) > 1
      bscalarflag = 1;
    end
  else
    if subs1.func.size(1)> 1 || subs1.func.size(2) > 1 ||...
        subs2.func.size(1)> 1 || subs2.func.size(2) > 1
      bscalarflag = 1;
    end
  end
end
if bscalarflag
  bMrow = bOver.func.size(1);
  bNcol = bOver.func.size(2);
  b.func.size = [bMrow,bNcol];
  % bOver has all derivatives such that they are the RepMatted ones, but
  % the repmatted derivatives have not yet been printed to file, before we
  % send b and bOver to cadaPrintReMap (to make sure they have same
  % OverMap), we need to RepMat the derivatives in b.
  for Vcount = 1:NUMvod
    if ~isempty(b.deriv(Vcount).name)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      brows = b.deriv(Vcount).nzlocs(:,1); bcols = b.deriv(Vcount).nzlocs(:,2);
      db = sparse(brows,bcols,1:length(brows),1,nv);
      db = repmat(db,bMrow*bNcol,1);
      [brows,bcols,binds] = find(db);
      if size(brows,2) > 1; brows = brows.'; bcols = bcols.'; end
      Dind1 = cadaindprint(binds);
      TDname = sprintf(['cada',NDstr,'td%1.0d'],Vcount);
      fprintf(fid,[indent,TDname,' = ',b.deriv(Vcount).name,'(',Dind1,');\n']);
      b.deriv(Vcount).name = TDname;
      b.deriv(Vcount).nzlocs = [brows,bcols];
    end
  end
end
% Make Sure that B is the same as the OverMapped B.
if isempty(b.func.value)
  b = cadaPrintReMap(b,bOver,b.id);
end

if isempty(b.func.name)
  b.func.name = btemp.func.name;
end

if isinf(y.func.size(1))
  yvec = 1;
elseif isinf(y.func.size(2))
  yvec = 2;
else
  yvec = 0;
end

CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
if DPFLAG == 1
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      y.deriv(Vcount).name = cadadername(funcstr,Vcount);
      if DPFLAG && ~isempty(b.deriv(Vcount).nzlocs)
        BindName = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES{Vcount,1};
        BindFlags = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES{Vcount,3};
        % -------Assign B derivatives---------%
        if ~isempty(BindName)
          if yvec
            % Vectorized
            if BindFlags(1)
              % Indices change on this loop
              if BindFlags(2)
                % Indices are Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,'(:,logical(',BindName,'(:,',CountName,'))) = ',...
                  b.deriv(Vcount).name,'(:,nonzeros(',BindName,'(:,',CountName,')));\n']);
              else
                % Indices are non-Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,' = ',...
                  b.deriv(Vcount).name,'(:,',BindName,'(:,',CountName,'));\n']);
              end
            else
              % Indices do not change on this loop
              if BindFlags(2)
                % Indices are Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,'(:,logical(',BindName,')) = ',...
                  b.deriv(Vcount).name,'(:,nonzeros(',BindName,'));\n']);
              else
                % Indices are non-Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,' = ',...
                  b.deriv(Vcount).name,'(:,',BindName,');\n']);
              end
            end
          else
            % Non-Vectorized
            if BindFlags(1)
              % Indices change on this loop
              if BindFlags(2)
                % Indices are Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,'(logical(',BindName,'(:,',CountName,'))) = ',...
                  b.deriv(Vcount).name,'(nonzeros(',BindName,'(:,',CountName,')));\n']);
              else
                % Indices are non-Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,' = ',...
                  b.deriv(Vcount).name,'(',BindName,'(:,',CountName,'));\n']);
              end
            else
              % Indices do not change on this loop
              if BindFlags(2)
                % Indices are Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,'(logical(',BindName,')) = ',...
                  b.deriv(Vcount).name,'(nonzeros(',BindName,'));\n']);
              else
                % Indices are non-Sparse
                fprintf(fid,[indent,y.deriv(Vcount).name,' = ',...
                  b.deriv(Vcount).name,'(',BindName,');\n']);
              end
            end
          end
        end
        
      end
      % -------Assign Zero derivatives---------%
      ZindName = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES{Vcount,4};
      ZindFlags = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).INDICES{Vcount,6};
      if ~isempty(ZindName)
        if ZindFlags(2)
          if yvec
            if ZindFlags(1)
              % Indices change on this loop
              fprintf(fid,[indent,y.deriv(Vcount).name,'(:,logical(',ZindName,'(:,',CountName,'))) = 0;\n']);
            else
              % Indices do not change on this loop
              fprintf(fid,[indent,y.deriv(Vcount).name,'(:,logical(',ZindName,')) = 0;\n']);
            end
          else
            if ZindFlags(1)
              % Indices change on this loop
              fprintf(fid,[indent,y.deriv(Vcount).name,'(logical(',ZindName,'(:,',CountName,'))) = 0;\n']);
            else
              % Indices do not change on this loop
              fprintf(fid,[indent,y.deriv(Vcount).name,'(logical(',ZindName,')) = 0;\n']);
            end
          end
        else
          % Indices are Non-Sparse => all zeros
          fprintf(fid,[indent,y.deriv(Vcount).name,'(:) = 0;\n']);
        end
      end
    end
  end
end


%-------------------Function Printing---------------%
if PFLAG == 1
  if ~isempty(x.func.name) && ~strcmp(funcstr,x.func.name)
    fprintf(fid,[indent,funcstr,' = ',x.func.name,';\n']);
  end
  if isempty(subs2) && ~isempty(ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{1,1})
    % Linear Assignment into a Matrix whose row size is changing - need to
    % fix the Assignment indices.
    RowIndName = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{1,1};
    RowDepFlag = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{1,3}(1);
    TFind1 = ['cada',NDstr,'tfind1']; TFind2 = ['cada',NDstr,'tfind2']; TFind3 = ['cada',NDstr,'tfind3'];
    fprintf(fid,[indent,TFind2,' = zeros(%1.0d,%1.0d);\n'],x.func.size(1),x.func.size(2));
    fprintf(fid,[indent,TFind2,'(:) = 1:%1.0d;\n'],x.func.size(1)*x.func.size(2));
    if RowDepFlag
      fprintf(fid,[indent,TFind3,' = ',TFind2,'(1:',RowIndName,'(',CountName,'),:);\n']);
    else
      fprintf(fid,[indent,TFind3,' = ',TFind2,'(1:',RowIndName,',:);\n']);
    end
    fprintf(fid,[indent,TFind1,' = ',TFind2,'(',subs1.func.name,');\n']);
    subs1.func.name = TFind1;
    % Set Logic Flag so that next section doesnt mess with the referencing
    % we just did.
    ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(1) = 2;
  end

  RowIndName = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{2,1}; % Size of B rows
  if ~isempty(RowIndName)
    RowDepFlag = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{2,3}(1);
    if RowDepFlag
      RowIndName = [RowIndName,'(',CountName,')'];
    end
  end
  ColIndName = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{3,1}; % Size of B cols
  if ~isempty(ColIndName)
    ColDepFlag = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).SIZES{3,3}(1);
    if ColDepFlag
      ColIndName = [ColIndName,'(',CountName,')'];
    end
  end
  RowLogicFlag = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(1);
  ColLogicFlag = ADIGATORFORDATA(INNERLOC).SUBSASGN(ASGNCOUNT).FLAGS(2);

  % ---Get RHS of what we are going to print.--- %
  if (isempty(RowIndName) && isempty(ColIndName)) || bscalarflag
    % B is scalar, or doesn't change sizes.
    RHSstr = [b.func.name,';\n'];
  elseif isempty(RowIndName)
    % B changing column sizes
    RHSstr = [b.func.name,'(:,1:',ColIndName,');\n'];
  elseif isempty(ColIndName)
    % B changing row sizes
    RHSstr = [b.func.name,'(1:',RowIndName,',:);\n'];
  else
    % B changing row and column sizes
    RHSstr = [b.func.name,'(1:',RowIndName,',1:',ColIndName,');\n'];
  end
  
  % ---Get LHS of what we are going to print.--- %
  if isempty(subs2)
    % Linear Indexing.
    if (isempty(RowIndName) && isempty(ColIndName)) || RowLogicFlag==2 || strcmp(subs1.func.name,':');
      % Either B doesnt change size, or the index is logical or ':'
      LHSstr = [funcstr,'(',subs1.func.name,')'];
    elseif isempty(RowIndName)
      % B changes column size
      if b.func.size(1) > 1
        LHSstr = sprintf([funcstr,'(',subs1.func.name,'(1:',ColIndName,'*%1.0f))'],b.func.size(1));
      else
        LHSstr = [funcstr,'(',subs1.func.name,'(1:',ColIndName,'))'];
      end
    elseif isempty(ColIndName)
      % B changes row size
      if b.func.size(2) > 1
        LHSstr = sprintf([funcstr,'(',subs1.func.name,'(1:',RowIndName,'*%1.0f))'],b.func.size(2));
      else
        LHSstr = [funcstr,'(',subs1.func.name,'(1:',RowIndName,'))'];
      end
    else
      % B changes row and column size
      if subs1.func.size(1) == 1 || subs1.func.size(2) == 1
        LHSstr = [funcstr,'(',subs1.func.name,'(1:',RowIndName,'*',ColIndName,'))'];
      else
        LHSstr = [funcstr,'(',subs1.func.name,'(1:',RowIndName,',','1:',ColIndName,'))'];
      end
    end
  else
    % Subs Referencing
    if isempty(RowIndName) || RowLogicFlag == 2
      Ref1 = subs1.func.name;
    elseif strcmp(subs1.func.name,':')
      Ref1 = ['1:',RowIndName];
    else
      Ref1 = [subs1.func.name,'(1:',RowIndName,')'];
    end
    if isempty(ColIndName) || ColLogicFlag == 2
      Ref2 = subs2.func.name;
    elseif strcmp(subs2.func.name,':')
      Ref2 = ['1:',ColIndName];
    else
      Ref2 = [subs2.func.name,'(1:',ColIndName,')'];
    end
    LHSstr = [funcstr,'(',Ref1,',',Ref2,')'];
  end
  fprintf(fid,[indent,LHSstr,' = ',RHSstr]);
end

if ~isempty(subs2)
  ADIGATOR.VARINFO.LASTOCC(subs2.id,1) = ADIGATOR.VARINFO.COUNT;
end
ADIGATOR.VARINFO.LASTOCC([subs1.id y.id x.id b.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(y,'cada')
  y = class(y,'cada');
end
return
end

function lRefStr = getLogicalDerivRef...
  (logicflag1,logicflag2,s,subs1,subs2,nsubs,FMrow,FNcol,nv)
global ADIGATOR
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
fid = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
lRefStr = ['cada',NDstr,'tf1'];
if logicflag1 && logicflag2 && ~isinf(FMrow) && ~isinf(FNcol)
  % Both Logical Unknowns - need to convert them to a single
  % reference.
  dim1 = numel(s.subs{1});
  dim2 = numel(s.subs{2});
  refind = 1:dim1*dim2;
  refind = reshape(refind,dim1,dim2);
  refind = refind(s.subs{1},s.subs{2});
  [refind1, refind2] = ind2sub(refind(:),dim1,dim2);
  refind1str = cadaindprint(refind1);
  refind2str = cadaindprint(refind2);
  fprintf(fid,[indent,lRefStr,' = ',subs1.func.name,'(',refind1str,') & ',...
    subs2.func.name,'(',refind2str,');\n']);
elseif logicflag1
  % First ref logical
  if nnz(s.subs{1}) == nsubs
    % Assignment locations fully defined by logical index
    if numel(s.subs{1}) == nsubs
      % No known zeros in logical index
      lRefStr = subs1.func.name;
    else
      % Have some known zeros in logical index.
      knownInd = cadaindprint(s.subs{1});
      fprintf(fid,[indent,lRefStr,' = ',subs1.func.name,'(',knownInd,');\n']);
      if isinf(FNcol)
        % Second dimension is vectorized - do this a little differently
        fprintf(fid,[indent,lRefStr,' = repmat(',lRefStr,'(:).'',[1 %1.0f]);\n'],nv);
      end
    end
  else
    % Second assignment index is vector
    if islogical(s.subs{2})
      dim = nnz(s.subs{2});
    elseif isnumeric(s.subs{2})
      dim = numel(s.subs{2});
    else
      dim = FNcol;
    end
    if nnz(s.subs{1}) < numel(s.subs{1})
      knownInd = cadaindprint(s.subs{1});
      fprintf(fid,[indent,lRefStr,' = ',subs1.func.name,'(',knownInd,');\n']);
      fprintf(fid,[indent,lRefStr,' = repmat(',lRefStr,',[1 %1.0f]);\n'],dim);
    else
      fprintf(fid,[indent,lRefStr,' = repmat(',subs1.func.name,'(:),[1 %1.0f]);\n'],dim);
    end
  end
else
  % Second ref logical
  if nnz(s.subs{2}) == nsubs
    % Assignment locations fully defined by logical index
    if numel(s.subs{2}) == nsubs
      % No known zeros in logical index
      lRefStr = subs2.func.name;
    else
      % Have some known zeros in logical index.
      knownInd = cadaindprint(s.subs{2});
      fprintf(fid,[indent,lRefStr,' = ',subs2.func.name,'(',knownInd,');\n']);
      if isinf(FMrow)
        % First Dimension is vectorized
        fprintf(fid,[indent,lRefStr,' = repmat(',lRefStr,'(:).'',[1 %1.0f]);\n'],nv);
      end
    end
  else
    % First assignment index is vector
    if islogical(s.subs{1})
      dim = nnz(s.subs{1});
    elseif isnumeric(s.subs{1})
      dim = numel(s.subs{1});
    else
      dim = FMrow;
    end
    if nnz(s.subs{2}) < numel(s.subs{2})
      knownInd = cadaindprint(s.subs{2});
      fprintf(fid,[indent,lRefStr,' = ',subs2.func.name,'(',knownInd,');\n']);
      fprintf(fid,[indent,lRefStr,' = repmat(',lRefStr,'(:).'',[%1.0f 1]);\n'],dim);
    else
      fprintf(fid,[indent,lRefStr,' = repmat(',subs2.func.name,'(:).'',[%1.0f 1]);\n'],dim);
    end
  end
end
return
end