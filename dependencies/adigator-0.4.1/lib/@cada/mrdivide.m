function z = mrdivide(x,y)
% CADA overloaded version of function MRDIVIDE
% If x and y are matrices, then z = x/y, and this function calls
% cadamtimesderiv one or more times to compute the derivative of z. If x or
% y is a scalar, then z = x./y is called instead.
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
  z = x./y;
  return
elseif xNcol == yNcol
  FMrow = xMrow;
  FNcol = yMrow;
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
  z.func.value = x.func.value/y.func.value;
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
    ztemp = xtemp/ytemp;
    z.func.zerolocs = find(~ztemp(:));
    if length(z.func.zerolocs) == FMrow*FNcol
      z.func.zerolocs = [];
      z.func.value    = zeros(FMrow,FNcol);
    end
    ztemp = abs(ztemp);
  end
end

TF3 = ['cada',NDstr,'tf3'];
funcprint = 0;
% ------------------------Build Derivative Properties----------------------
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod;
  % For the Derivative of Z = X/Y, we have the rule:
  % dZ = dX/Y - [(X/Y)*dY]/Y = [dX - (X/Y)*dY]/Y
  dxinds = y.deriv(Vcount).nzlocs;
  if yMrow ~= yNcol && ~isempty(dxinds)
    error('mrdivide derivatives for non-square matrix not written yet');
  end
  if ~isempty(y.deriv(Vcount).nzlocs)
    if DPFLAG && ~funcprint
      fprintf(fid,[indent,TF3,' = -',x.func.name,'/',y.func.name,';\n']);
      funcprint = 1;
    end
    % --------------------------- -(X/Y)*dY ----------------------------- %
    % Let A have the function properties of -(X/Y) and no derivatives
    A.func      = z.func;
    A.func.name = TF3;
    Atemp       = ztemp;
    A.deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    TD3 = ['cada',NDstr,'td3'];
    % Let B = A*Y, dB = A*dY = -(X/Y)*dY
    Bnzlocs = cadamtimesderiv(A,x,Atemp,xtemp,Vcount,TD3,DPFLAG,'mtimes');
    % dim(A) = dim(z), so, dim(A*Y) = [zMrow yNcol] = dim(X)
    B.func = y.func; B.func.name = []; Btemp = ytemp; %Btemp never really gets used
    B.deriv = A.deriv;
    B.deriv(Vcount).name   = TD3;
    B.deriv(Vcount).nzlocs = Bnzlocs;
  else
    Bnzlocs = [];
  end
  
  if ~isempty(Bnzlocs) && ~isempty(dxinds)
    % ------------------------ dX - (X/Y)*dY ---------------------------- %
    % need to union derivative the two derivatives prior to performing
    % the mldivide
    [xBnzlocs,dBinds,dXinds,Bindflag,Xindflag] = cadaunion(Bnzlocs,...
      dxinds,xMrow*xNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
    TD4 = ['cada',NDstr,'td4'];
    if DPFLAG
      if Bindflag && Xindflag
        fprintf(fid,[indent,TD4,' = ',TD3,' + ',x.deriv(Vcount).name,';\n']);
      else
        if Bindflag
          fprintf(fid,[indent,TD4,' = ',TD3,';\n']);
        else
          fprintf(fid,[indent,TD4,' = zeros(%1.0f,1);\n'],size(xBnzlocs,1));
          Dind1 = cadaindprint(dBinds);
          fprintf(fid,[indent,TD4,'(',Dind1,') = ',TD3,';\n']);
        end
        if Xindflag
          fprintf(fid,[indent,TD4,' = ',TD4,' + ',x.deriv(Vcount).name,';\n']);
        else
          Dind1 = cadaindprint(dXinds);
          fprintf(fid,[indent,TD4,'(',Dind1,') = ',TD4,'(',Dind1,') + ',x.deriv(Vcount).name,';\n']);
        end
      end
    end
    % Let C = X + B, dC = dX + dB = dX -(X/Y)*dY
    C.func = x.func; C.func.name = []; Ctemp = xtemp;
    C.deriv = A.deriv;
    C.deriv(Vcount).name   = TD4;
    C.deriv(Vcount).nzlocs = xBnzlocs;
    
    % ----------------------- [dX - (X/Y)*dY]/Y ------------------------- %
    % let D have the function properties of Y with no derivatives
    D.func  = y.func; Dtemp = ytemp;
    D.deriv = A.deriv;
    % Z = C/D, dZ = dC/D = [dX - (X/Y)*dY]/Y
    derivstr = cadadername(funcstr,Vcount);
    nzlocs = cadamtimesderiv(C,D,Ctemp,Dtemp,Vcount,derivstr,DPFLAG,'mrdivide');
    if ~isempty(nzlocs)
      z.deriv(Vcount).name   = derivstr;
      z.deriv(Vcount).nzlocs = nzlocs;
    end
    
  elseif ~isempty(Bnzlocs)
    % ---------------------------- -[(X/Y)*dY]/Y ------------------------ %
    % let D have the function properties of Y with no derivatives
    D.func  = y.func; Dtemp = ytemp;
    D.deriv = A.deriv;
    % Z = B/D, so dZ = dB/D = -[(X/Y)*dY]/Y
    derivstr = cadadername(funcstr,Vcount);
    nzlocs = cadamtimesderiv(B,D,Btemp,Dtemp,Vcount,derivstr,DPFLAG,'mrdivide');
    if ~isempty(nzlocs)
      z.deriv(Vcount).name   = derivstr;
      z.deriv(Vcount).nzlocs = nzlocs;
    end
  
  elseif ~isempty(dxinds)
    % ----------------------------- dX/Y -------------------------------- %
    if ~isempty(y.deriv(Vcount).nzlocs)
      % for some reason y had derivs but B didnt, some odd zero
      % cancellation
      y.deriv(Vcount).nzlocs = []; y.deriv(Vcount).name = [];
    end
    derivstr = cadadername(funcstr,Vcount);
    nzlocs = cadamtimesderiv(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,'mrdivide');
    if ~isempty(nzlocs)
      z.deriv(Vcount).name   = derivstr;
      z.deriv(Vcount).nzlocs = nzlocs;
    end
    
  end
    
end

% --------------------------Function Printing --------------------------- %

if PFLAG == 1 && ~funcprint
  fprintf(fid,[indent,funcstr,' = ',x.func.name,'/',y.func.name,';\n']);
elseif funcprint
  fprintf(fid,[indent,funcstr,' = -',TF3,';\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;