function z = mldivide(x,y)
% CADA overloaded version of function MLDIVIDE
% If x and y are matrices, then z = x\y, and this function calls
% cadamtimesderiv one or more times to compute the derivative of z. If x or
% y is a scalar, then z = y./x is called instead.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  z = cadaEmptyEval(x,y);
  return
end
PFLAG  = ADIGATOR.PRINT.FLAG;
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% ----------------------------Parse Inputs------------------------------- %
if isa(x,'cada') && isa(y,'cada')
  % Both Inputs are Symbolic
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  yMrow = y.func.size(1); yNcol = y.func.size(2);
elseif isa(x,'cada')
  % y is numeric input
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  [yMrow,yNcol] = size(y);
  ytemp.id = [];
  ytemp.func = struct('name',[],'size',[yMrow yNcol],'zerolocs',[],...
    'value',y);
  if PFLAG
    if yMrow*yNcol == 1
      ytemp.func.name = num2str(y,16);
    else
      ytemp.func.name = cadamatprint(y);
    end
  end
  ytemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  y = ytemp;
  y = class(y,'cada');
else
  % x is numeric input
  yMrow = y.func.size(1); yNcol = y.func.size(2);
  [xMrow,xNcol] = size(x);
  xtemp.id = [];
  xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
    'value',x);
  if PFLAG
    if xMrow*xNcol == 1
      xtemp.func.name = num2str(x,16);
    else
      xtemp.func.name = cadamatprint(x);
    end
  end
  xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  x = xtemp;
  x = class(x,'cada');
end

% ----------------------------Function Sizing------------------------------
if (xMrow == 1 && xNcol == 1) || (yMrow == 1 && yNcol == 1)
  z = y./x;
  return
elseif xMrow == yMrow
  FMrow = xNcol;
  FNcol = yNcol;
  if isinf(FMrow) || isinf(FNcol)
    error('Cannot mrdivide if matrix has vectorized dimension')
  end
else
  error('Inputs are not of compatible sizes');
end

%-------------------------------------------------------------------------%
%                      Build Function Properties                          %
%-------------------------------------------------------------------------%
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);
if ~isempty(x.func.value) && ~isempty(y.func.value)
  % z is numeric
  z.func.value = x.func.value\y.func.value;
else
  spflag = 0;
  if ~isempty(x.func.value)
    xtemp  = x.func.value;
    spflag = 1;
  elseif ~isempty(x.func.zerolocs)
    xtemp = rand(xMrow,xNcol);
    xtemp(x.func.zerolocs) = 0;
    spflag = 1;
  else
    xtemp = rand(xMrow,xNcol);
  end
  if ~isempty(y.func.value)
    ytemp  = y.func.value;
    spflag = 1;
  elseif ~isempty(y.func.zerolocs)
    ytemp = rand(yMrow,yNcol);
    ytemp(y.func.zerolocs) = 0;
    spflag = 1;
  else
    ytemp = rand(yMrow,yNcol);
  end
  if spflag == 1
    ztemp = xtemp\ytemp;
    z.func.zerolocs = find(~ztemp(:));
    if length(z.func.zerolocs) == FMrow*FNcol
      z.func.zerolocs = [];
      z.func.value    = zeros(FMrow,FNcol);
    end
    ztemp = abs(ztemp);
  else
    ztemp = abs(xtemp\ytemp);
  end
end

TF3 = ['cada',NDstr,'tf3'];
TF4 = ['cada',NDstr,'tf4'];
funcprint = 0;
% ------------------------Build Derivative Properties----------------------
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod;
  if xMrow == xNcol
    % ----------------------- SQUARE SYSTEM ----------------------------- %
    % for this system we have that X*Z - Y = 0, thus dX*Z + X*dZ - dY = 0;
    % so, dZ = -X\(dX*Z) + X\dY
    % For the Derivative of Z = X\Y, we have the rule:
    % dZ = -X\[dX*(X\Y)] + X\dY = X\[-dX*(X\Y) + dY]
    dyinds = y.deriv(Vcount).nzlocs;
    
    if ~isempty(x.deriv(Vcount).nzlocs)
      if DPFLAG && ~funcprint
        fprintf(fid,[indent,TF3,' = -',x.func.name,'\\',y.func.name,';\n']);
        funcprint = 1;
      end
      % --------------------------- -dX*(X\Y) ----------------------------- %
      % Let A have the function properties of (X\Y) and no derivatives
      % A = (X\Y); dA = 0;
      A.func      = z.func;
      A.func.name = TF3;
      Atemp       = ztemp;
      A.deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      TD3 = ['cada',NDstr,'td3'];
      % Let B = X*A, dB = dX*(X\Y)
      Bnzlocs = cadamtimesderiv(x,A,xtemp,Atemp,Vcount,TD3,DPFLAG,'mtimes');
      % dim(A) = dim(z), so, dim(X*A) = [xMrow FNcol] = dim(Y)
      B.func = y.func; B.func.name = []; Btemp = ytemp; %Btemp never really gets used
      B.deriv = A.deriv;
      B.deriv(Vcount).name   = TD3;
      B.deriv(Vcount).nzlocs = Bnzlocs;
    else
      Bnzlocs = [];
    end
    
    if ~isempty(Bnzlocs) && ~isempty(dyinds)
      % ---------------------- -dX*(X\Y) + dY --------------------------- %
      % need to union derivative the two derivatives prior to performing
      % the mldivide
      [yBnzlocs,dBinds,dYinds,Bindflag,Yindflag] = cadaunion(Bnzlocs,...
        dyinds,yMrow*yNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
      TD4 = ['cada',NDstr,'td4'];
      if DPFLAG
        if Bindflag && Yindflag
          fprintf(fid,[indent,TD4,' = ',TD3,' + ',y.deriv(Vcount).name,';\n']);
        else
          if Bindflag
            fprintf(fid,[indent,TD4,' = ',TD3,';\n']);
          else
            fprintf(fid,[indent,TD4,' = zeros(%1.0f,1);\n'],size(yBnzlocs,1));
            Dind1 = cadaindprint(dBinds);
            fprintf(fid,[indent,TD4,'(',Dind1,') = ',TD3,';\n']);
          end
          if Yindflag
            fprintf(fid,[indent,TD4,' = ',TD4,' + ',y.deriv(Vcount).name,';\n']);
          else
            Dind1 = cadaindprint(dYinds);
            fprintf(fid,[indent,TD4,'(',Dind1,') = ',TD4,'(',Dind1,') + ',y.deriv(Vcount).name,';\n']);
          end
        end
      end
      % Let C = B + Y, dC = -dX*(X\Y) + dY
      C.func = y.func; C.func.name = []; Ctemp = ytemp;
      C.deriv = A.deriv;
      C.deriv(Vcount).name   = TD4;
      C.deriv(Vcount).nzlocs = yBnzlocs;
      
      % -------------------- X\[-dX*(X\Y) + dY] ------------------------- %
      % let D have the function properties of X with no derivatives
      D.func  = x.func; Dtemp = xtemp;
      D.deriv = A.deriv;
      % Z = D\C, dZ = D\dC = X\[-dX*(X\Y) + dY]
      derivstr = cadadername(funcstr,Vcount);
      nzlocs = cadamtimesderiv(D,C,Dtemp,Ctemp,Vcount,derivstr,DPFLAG,'mldivide');
      if ~isempty(nzlocs)
        z.deriv(Vcount).name   = derivstr;
        z.deriv(Vcount).nzlocs = nzlocs;
      end
      
    elseif ~isempty(Bnzlocs)
      % ---------------------------- -X\[dX*(X\Y)] ------------------------ %
      % let D have the function properties of X with no derivatives
      D.func  = x.func; Dtemp = xtemp;
      D.deriv = A.deriv;
      % Z = D\B, dZ = D\dB = -X\[dX*(X\Y)]
      derivstr = cadadername(funcstr,Vcount);
      nzlocs = cadamtimesderiv(D,B,Dtemp,Btemp,Vcount,derivstr,DPFLAG,'mldivide');
      if ~isempty(nzlocs)
        z.deriv(Vcount).name   = derivstr;
        z.deriv(Vcount).nzlocs = nzlocs;
      end
      
    elseif ~isempty(dyinds)
      % ----------------------------- X\dY -------------------------------- %
      if ~isempty(x.deriv(Vcount).nzlocs)
        % for some reason x had derivs but B didnt, some odd zero
        % cancellation
        x.deriv(Vcount).nzlocs = []; x.deriv(Vcount).name = [];
      end
      derivstr = cadadername(funcstr,Vcount);
      nzlocs = cadamtimesderiv(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,'mldivide');
      if ~isempty(nzlocs)
        z.deriv(Vcount).name   = derivstr;
        z.deriv(Vcount).nzlocs = nzlocs;
      end
      
    end
  elseif xMrow > xNcol
    % ------------------- OVERDETERMINED SYSTEM ------------------------- %
    % For this system we have that X'*(X*Z-Y) = 0, thus 
    % (X*Z-Y)'*dX + X'*(dX*Z + X*dZ - dY) = 0;
    % then X'*X*dZ = - (X*Z-Y)'*dX - X'*dX*Z + X'*dY
    % where (X'*X) must have rank xNcol if X has rank xNcol so it is
    % invertible
    % Thus,
    % dZ = -(X'*X)\[(X*Z-Y)'*dX + X'*(dX*Z - *dY)]
    if ~isempty(x.deriv(Vcount).nzlocs)
      if DPFLAG && ~funcprint
        fprintf(fid,[indent,TF3,' = ',x.func.name,'\\',y.func.name,';\n']);
        funcprint = 2;
      end
      %dZ = -(X'*X)\[((X*Z-Y)'*dX)' + X'*dX*Z]
      % --------------------------- X'*dX*Z ----------------------------- %
      
      % First lets define A to have the function properties of Z with no
      % derivatives
      % A = (X\Y), dA = 0;
      A.func      = z.func;
      A.func.name = TF3;
      Atemp       = ztemp;
      A.deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      % want to print out the calculations for dX*Z, so let 
      % B = X*A, dB = dX*A = dX*Z;
      TD3 = ['cada',NDstr,'td3'];
      Bnzlocs = cadamtimesderiv(x,A,xtemp,Atemp,Vcount,TD3,DPFLAG,'mtimes');
      
      % dim(A) = dim(z), so, dim(X*A) = [xMrow FNcol] = dim(Y)
      if ~isempty(Bnzlocs)
        B.func = y.func; B.func.name = []; Btemp = ytemp; %Btemp never really gets used
        B.deriv = A.deriv;
        B.deriv(Vcount).name   = TD3;
        B.deriv(Vcount).nzlocs = Bnzlocs;
        % Now we define C to have the function properties of X' with no
        % derivatives, so C = X', dC = 0;
        if DPFLAG
          fprintf(fid,[indent,TF4,' = ',x.func.name,'.'';\n']);
        end
        C.func = x.func;
        C.func.size = [xNcol xMrow];
        C.func.name = TF4;
        Ctemp  = xtemp.';
        C.deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        % We want to print out the calculations for X'*dX*Z, so let
        % D = C*B, dC = C*dB = X'*dX*Z
        TD4 = ['cada',NDstr,'td4'];
        Dnzlocs = cadamtimesderiv(C,B,Ctemp,Btemp,Vcount,TD4,DPFLAG,'mtimes');
      else
        Dnzlocs = [];
      end
      % -------------------------- (X*Z-Y)'*dX -------------------------- %
      % Let E have the function properties of (X*Z-Y)' and no derivatives,
      % E = (X*Z-Y)'; dE = 0;
      if DPFLAG
        fprintf(fid,[indent,TF4,' = (',x.func.name,'*',TF3,' - ',y.func.name,').'';\n']);
      end
      E.func      = y.func;
      E.func.size = [yNcol yMrow];
      E.func.name = TF4;
      Etemp       = (xtemp*ztemp-ytemp).';
      E.deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      % We want to print out the calculations for (X*Z-Y)'*dX, so let 
      % F = E*X, thus dF = E*dX = (X*Z-Y)'*dX
      TD5 = ['cada',NDstr,'td5'];
      Fnzlocs = cadamtimesderiv(E,x,Etemp,xtemp,Vcount,TD5,DPFLAG,'mtimes');
      
      if ~isempty(Dnzlocs) && ~isempty(Fnzlocs) && isempty(y.deriv(Vcount).nzlocs)
        % ---------------- ((X*Z-Y)'*dX)' + X'*dX*Z --------------------- %
        %                       dF' + dD
        D.func = x.func;
        D.func.name = [];
        D.func.size = [C.func.size(1) B.func.size(2)];
        D.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        D.deriv(Vcount).nzlocs = Dnzlocs;
        D.deriv(Vcount).name   = TD4;
        
        F.func = x.func;
        F.func.name = [];
        F.func.size = [E.func.size(1) x.func.size(2)];
        F.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        F.deriv(Vcount).nzlocs = Fnzlocs;
        F.deriv(Vcount).name   = TD5;
        
        % Need to transpose derivative of if F isnt a vector
        if F.func.size(1) > 1 && F.func.size(2) > 1
          fTranMap = zeros(F.func.size);
          fTranMap(:) = 1:F.func.size(1)*F.func.size(2);
          fTranMap = fTranMap.';
          nv = ADIGATOR.VAROFDIFF(Vcount).usize;
          dF = sparse(Fnzlocs(:,1),Fnzlocs(:,2),1:size(Fnzlocs,1),...
            F.func.size(1)*F.func.size(2),nv);
          dFtran = dF(fTranMap(:),:);
          [ftrows, ftcols,ftlocs] = find(dFtran);
          if size(ftrows,2) > 1; ftrows = ftrows.'; ftcols = ftcols.'; end
          F.deriv(Vcount).nzlocs = [ftrows ftcols];
          if DPFLAG
            TD5inds = cadaindprint(ftlocs(:));
            fprintf(fid,[indent,TD5,' = ',TD5,'(',TD5inds,');\n']);
          end
        end
        F.func.size = D.func.size;
        
        % need to union derivative the D and F prior to performing
        % the mldivide
        % let G = D+F, G - TD3, D - TD4, F - TD5
        [Gnzlocs,dDinds,dFinds,Dindflag,Findflag] = cadaunion(Dnzlocs,...
          F.deriv(Vcount).nzlocs,FMrow*FNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
        if DPFLAG
          if Dindflag && Findflag
            fprintf(fid,[indent,TD3,' = ',TD4,' + ',TD5,';\n']);
          else
            if Dindflag
              fprintf(fid,[indent,TD3,' = ',TD4,';\n']);
            else
              fprintf(fid,[indent,TD3,' = zeros(%1.0f,1);\n'],size(Gnzlocs,1));
              Dind1 = cadaindprint(dDinds);
              fprintf(fid,[indent,TD3,'(',Dind1,') = ',TD4,';\n']);
            end
            if Findflag
              fprintf(fid,[indent,TD3,' = ',TD3,' + ',TD5,';\n']);
            else
              Dind1 = cadaindprint(dFinds);
              fprintf(fid,[indent,TD3,'(',Dind1,') = ',TD3,'(',Dind1,') + ',TD5,';\n']);
            end
          end
        end
        % Let G = D+F, dG = (X*Z-Y)'*dX)' + X'*dX*Z
        G.func = D.func;
        Gtemp = ztemp;
        G.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        G.deriv(Vcount).nzlocs = Gnzlocs;
        G.deriv(Vcount).name   = TD3;
        % Let H = -(X'*X), dH = 0;
        if DPFLAG
          fprintf(fid,[indent,TF4,' = -(',x.func.name,'.''*',x.func.name,');\n']);
        end
        H.func = x.func;
        H.func.size = [xNcol xNcol];
        H.func.name = TF4;
        H.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        Htemp = xtemp.'*xtemp;
        
        % now we compute z = H\G, dz = H\dG = -(X'*X)\[((X*Z-Y)'*dX)' + X'*dX*Z]
        derivstr = cadadername(funcstr,Vcount);
        nzlocs = cadamtimesderiv(H,G,Htemp,Gtemp,Vcount,derivstr,DPFLAG,'mldivide');
        if ~isempty(nzlocs)
          z.deriv(Vcount).nzlocs = nzlocs;
          z.deriv(Vcount).name   = derivstr;
        end
      else
        error('not coded yet')
      end
    end
  elseif xMrow < xNcol
    error('Not coded yet')
  end
end

% --------------------------Function Printing --------------------------- %

if PFLAG == 1 && ~funcprint
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'\\',y.func.name,';\n']);
elseif funcprint == 1
  fprintf(fid,[indent,funcstr,' = -',TF3,';\n']);
elseif funcprint == 2
  fprintf(fid,[indent,funcstr,' = ',TF3,';\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;