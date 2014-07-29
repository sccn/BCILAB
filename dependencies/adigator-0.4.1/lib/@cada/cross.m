function z = cross(x,y,varargin)
% CADA cross product

global ADIGATOR

if ADIGATOR.EMPTYFLAG
  if nargin > 2
    z = cadaEmptyEval(x,y,varargin{1});
  else
    z = cadaEmptyEval(x,y);
  end
  return
end

NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Parse Inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
if isa(x,'cada') && isa(y,'cada')
  % Both Inputs are Symbolic
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  yMrow = y.func.size(1); yNcol = y.func.size(2);
  if isinf(xMrow); xvec = 1; elseif isinf(xNcol); xvec =2; else xvec = 0;end
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

%% ~~~~~~~~~~~~~~~~~~~~~~~~ Function Sizing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
dumbflag = 0;
if (xMrow == yMrow && xNcol == yNcol)
  FMrow = yMrow; FNcol = yNcol;
elseif xMrow ==1 && xNcol == 3 && yMrow == 3 && yNcol == 1 && nargin == 2
  dumbflag = 1;
  %yMrow = 1; yNcol = 3;
  FMrow = 1; FNcol = 3;
elseif xMrow ==3 && xNcol == 1 && yMrow == 1 && yNcol == 3 && nargin == 2
  dumbflag = 1;
  %xMrow = 1; xNcol = 3;
  FMrow = 1; FNcol = 3;
else
  error('Inputs must be same size')
end
if nargin > 2
  dim = varargin{1};
  if isa(dim,'cada')
    if ~isempty(dim.func.value)
      dim = dim.func.value;
    else
      error('Dimension to cross must be known numeric')
    end
  end
  if (dim == 1 && FMrow ~= 3) || (dim == 2 && FNcol ~= 3) || ~any(dim == [1 2])
    error('Invalid dimension input')
  end
elseif FMrow == 3
  dim = 1;
elseif FNcol == 3
  dim = 2;
else
  error('No input dimension is 3')
end

%% ~~~~~~~~~~~~~~~~~~~~ Build Function Properties ~~~~~~~~~~~~~~~~~~~~~~ %%
z.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
z.func = struct('name',funcstr,'size',[FMrow FNcol],'zerolocs',[],...
  'value',[]);

% Get Rid of vectorized dims
if isinf(xMrow); xMrow = 1; elseif isinf(xNcol); xNcol = 1; end
if isinf(yMrow); yMrow = 1; elseif isinf(yNcol); yNcol = 1; end
if isinf(FMrow); FMrow = 1; elseif isinf(FNcol); FNcol = 1; end

% ----------------------- Get I/J/K Reference Indices ------------------- %
refInds = reshape(1:FMrow*FNcol,[FMrow, FNcol]);
toughref = 0;
if dim == 1
  iref = refInds(1,:).';  
  jref = refInds(2,:).';
  kref = refInds(3,:).';
  if FNcol > 1
    toughref = 1;
    tranref = refInds.';
    tranref = tranref(:);
  end
else
  iref = refInds(:,1);
  jref = refInds(:,2);
  kref = refInds(:,3);
end
% Using this notation:
%
%   | xi |        | yi |    | zi |   | xj*yk - xk*yj |
%   | xj | cross  | yj |  = | zj | = | xk*yi - xi*yk |
%   | xk |        | yk |    | zk |   | xi*yj - xj*yi |

% --------------Function Numeric Values and Sparsity--------------------- %
if ~isempty(x.func.value) && ~isempty(y.func.value)
  % z is numeric
  z.func.value = cross(x.func.value,y.func.value,dim);
else
  spxflag = 0; spyflag = 0;
  if ~isempty(x.func.value)
    xtemp = logical(x.func.value); spxflag = 1;
    if dumbflag; xtemp = reshape(xtemp,FMrow,FNcol); end
  elseif ~isempty(x.func.zerolocs)
    xtemp = true(FMrow,FNcol); xtemp(x.func.zerolocs) = false; spxflag = 1;
  else
    xtemp = true(FMrow,FNcol);
  end
  if ~isempty(y.func.value)
    ytemp = logical(y.func.value); spyflag = 1;
    if dumbflag; ytemp = reshape(ytemp,FMrow,FNcol); end
  elseif ~isempty(y.func.zerolocs)
    ytemp = true(FMrow,FNcol); ytemp(y.func.zerolocs) = false; spyflag = 1;
  else
    ytemp = true(FMrow,FNcol);
  end
  if spxflag || spyflag
    xitemp = xtemp(iref); yitemp = ytemp(iref);
    xjtemp = xtemp(jref); yjtemp = ytemp(jref);
    xktemp = xtemp(kref); yktemp = ytemp(kref);
    ztemp  = false(FMrow,FNcol);
    ztemp(iref) = (xjtemp & yktemp) | (xktemp & yjtemp);
    ztemp(jref) = (xktemp & yitemp) | (xitemp & yktemp);
    ztemp(kref) = (xitemp & yjtemp) | (xjtemp & yitemp);
    z.func.zerolocs = find(~ztemp(:));
    if length(z.func.zerolocs) == FMrow*FNcol
      z.func.zerolocs = []; z.func.value = zeros(FMrow,FNcol);
    end
  end
  if spxflag
    fxLtemp = false(FMrow*FNcol,1);
    fxLtemp([kref;iref;jref]) = [xjtemp;xktemp;xitemp];
    fxRtemp = false(FMrow*FNcol,1);
    fxRtemp([jref;kref;iref]) = [xktemp;xitemp;xjtemp];
  end
  if spyflag
    fyLtemp = false(FMrow*FNcol);
    fyLtemp([jref;kref;iref]) = [yktemp;yitemp;yjtemp];
    fyRtemp = false(FMrow*FNcol,1);
    fyRtemp([kref;iref;jref]) = [yjtemp;yktemp;yitemp];
  end
end

%% ~~~~~~~~~~~~~~~~~~~~ Build Derivative Properties ~~~~~~~~~~~~~~~~~~~~ %%
% Define the following
%
%       | xj |        | xk |       | yk |       | yj |
% xL =  | xk |, xR  = | xi |, yL = | yi |, yR = | yk |
%       | xi |        | xj |       | yj |       | yi |
%     Z = xL.*yL - xR.*yR
% so dZ = dxL.*yL + xL.*dyL - dxR.*yR - xR.*dyR
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod;
  nv = ADIGATOR.VAROFDIFF(Vcount).usize;
  dxind = x.deriv(Vcount).nzlocs;
  dyind = y.deriv(Vcount).nzlocs;
  % ---------------Use Function Sparsity to Cancel Derivatives----------- %
  spdxLflag = 0; spdxRflag = 0; spdyLflag = 0; spdyRflag = 0;
  if ~isempty(dyind) && spxflag
    fxL = fxLtemp(dyind(:,1)); nzlocsyL = logical(fxL(:).*dyind(:,1));
    spdyLind = dyind(nzlocsyL,:);
    if size(spdyLind,1) < size(dyind,1); spdyLflag = 1; end
    fxR = fxRtemp(dyind(:,1)); nzlocsyR = logical(fxR(:).*dyind(:,1));
    spdyRind = dyind(nzlocsyR,:);
    if size(spdyRind,1) < size(dyind,1); spdyRflag = 1; end
  else
    spdyLind = dyind;
    spdyRind = dyind;
  end
  if ~isempty(dxind) && spyflag
    fyL = fyLtemp(dxind(:,1)); nzlocsxL = logical(fyL(:).*dxind(:,1));
    spdxLind = dxind(nzlocsxL,:);
    if size(spdxLind,1) < size(dxind,1); spdxLflag = 1; end
    fyR = fyRtemp(dxind(:,1)); nzlocsxR = logical(fyR(:).*dxind(:,1));
    spdxRind = dxind(nzlocsxR,:);
    if size(spdxRind,1) < size(dxind,1); spdxRflag = 1; end
  else
    spdxLind = dxind;
    spdxRind = dxind;
  end
  %% --------------------------- Build dz ------------------------------ %%
  dzflag = 0;
  if ~isempty(spdxLind) || ~isempty(spdxRind) || ...
      ~isempty(spdyLind) || ~isempty(spdyRind)
    % Initialize dz
    dz = sparse([],[],[],FMrow*FNcol,nv);
    dzflag = 1;
  end
  if ~isempty(spdxLind) % dxL
    nzxL = size(spdxLind,1);
    dxL  = sparse(spdxLind(:,1),spdxLind(:,2),1:nzxL,FMrow*FNcol,nv);
    dxL  = dxL([jref;kref;iref],:);
    if toughref; dxL = dxL(tranref,:); end
    dz = dz + dxL;
  end
  if ~isempty(spdxRind) % dxR
    nzxR = size(spdxRind,1);
    dxR  = sparse(spdxRind(:,1),spdxRind(:,2),1:nzxR,FMrow*FNcol,nv);
    dxR  = dxR([kref;iref;jref],:);
    if toughref; dxR = dxR(tranref,:); end
    dz = dz + dxR;
  end
  if ~isempty(spdyLind) % dyL
    nzyL = size(spdyLind,1);
    dyL  = sparse(spdyLind(:,1),spdyLind(:,2),1:nzyL,FMrow*FNcol,nv);
    dyL  = dyL([kref;iref;jref],:);
    if toughref; dyL = dyL(tranref,:); end
    dz = dz + dyL;
  end
  if ~isempty(spdyRind) % dyR
    nzyR = size(spdyRind,1);
    dyR  = sparse(spdyRind(:,1),spdyRind(:,2),1:nzyR,FMrow*FNcol,nv);
    dyR  = dyR([jref;kref;iref],:);
    if toughref; dyR = dyR(tranref,:); end
    dz = dz + dyR;
  end
  if dzflag
    % Get locations of dz
    derivstr = cadadername(funcstr,Vcount);
    z.deriv(Vcount).name = derivstr;
    [zrows,zcols] = find(dz);
    z.deriv(Vcount).nzlocs = [zrows,zcols];
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      nz = length(zrows);
      dz = sparse(zrows,zcols,1:nz,FMrow*FNcol,nv);
      dzinit = 0;
    end
  end
  %% --------------------- Print out dxL Calcs ------------------------- %%
  if DPFLAG && ~isempty(spdxLind) % dxL
    [~,dxLinds] = sort(nonzeros(dxL));
    % Get dxL string
    if spdxLflag
      dxLindin = (1:size(dxind,1)).'; 
      dxLindin = dxLindin(nzlocsxL);
      dxLref   = cadaindprint(dxLindin);
      TD2 = ['cada',NDstr,'td2'];
      if xvec
        fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,'(:',dxLref,');']);
      else
        fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,'(',dxLref,');']);
      end
      dxLstr = TD2;
    else
      dxLstr = x.deriv(Vcount).name;
    end
    % Get yL string
    yLref = [kref;iref;jref];
    if toughref; yLref = yLref(tranref); end
    yLref = nonzeros(diag(yLref)*logical(dxL));
    yLref = cadaindprint(yLref(dxLinds));
    if xvec == 1
      yLstr = [y.func.name,'(:,',yLref,')'];
    elseif xvec == 2
      yLstr = [y.func.name,'(',yLref,',:).'''];
    elseif yMrow == 1
      yLstr = [y.func.name,'(',yLref,').'''];
    else
      yLstr = [y.func.name,'(',yLref,')'];
    end
    % Get assignment to dz
    dzxLind  = full(dz(logical(dxL)));
    dzxLasgn = cadaindprint(dzxLind(dxLinds));
    if xvec
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nz);
      end
      dzxLstr  = [TD1,'(:,',dzxLasgn,')'];
    else
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nz);
      end
      dzxLstr  = [TD1,'(',dzxLasgn,')'];
    end
    % Print out the calculation
    if ~dzinit
      fprintf(fid,[indent,dzxLstr,' = ',dxLstr,'.*',yLstr,';\n']);
      dzinit = 1;
    else
      fprintf(fid,[indent,dzxLstr,' = ',dzxLstr,' + ',dxLstr,'.*',yLstr,';\n']);
    end
  end
  %% --------------------- Print out dxR Calcs ------------------------- %%
  if DPFLAG && ~isempty(spdxRind) % dxR
    [~,dxRinds] = sort(nonzeros(dxR));
    % Get dxR string
    if spdxRflag
      dxRindin = (1:size(dxind,1)).'; 
      dxRindin = dxRindin(nzlocsxR);
      dxRref   = cadaindprint(dxRindin);
      TD2 = ['cada',NDstr,'td2'];
      if xvec
        fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,'(:',dxRref,');']);
      else
        fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,'(',dxRref,');']);
      end
      dxRstr = TD2;
    else
      dxRstr = x.deriv(Vcount).name;
    end
    % Get yR string
    yRref = [jref;kref;iref];
    if toughref; yRref = yRref(tranref); end
    yRref = nonzeros(diag(yRref)*logical(dxR));
    yRref = cadaindprint(yRref(dxRinds));
    if xvec == 1
      yRstr = [y.func.name,'(:,',yRref,')'];
    elseif xvec == 2
      yRstr = [y.func.name,'(',yRref,',:).'''];
    elseif yMrow == 1
      yRstr = [y.func.name,'(',yRref,').'''];
    else
      yRstr = [y.func.name,'(',yRref,')'];
    end
    % Get assignment to dz
    dzxRind  = full(dz(logical(dxR)));
    dzxRasgn = cadaindprint(dzxRind(dxRinds));
    if xvec
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(size(',x.deriv(Vcount).name,',1),%1.0f);\n'],nz);
      end
      dzxRstr  = [TD1,'(:,',dzxRasgn,')'];
    else
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nz);
      end
      dzxRstr  = [TD1,'(',dzxRasgn,')'];
    end
    % Print out the calculation
    if ~dzinit
      fprintf(fid,[indent,dzxRstr,' = -',dxRstr,'.*',yRstr,';\n']);
      dzinit = 1;
    else
      fprintf(fid,[indent,dzxRstr,' = ',dzxRstr,' - ',dxRstr,'.*',yRstr,';\n']);
    end
  end
  %% --------------------- Print out dyL Calcs ------------------------- %%
  if DPFLAG && ~isempty(spdyLind) % dyL
    [~,dyLinds] = sort(nonzeros(dyL));
    % Get dyL string
    if spdyLflag
      dyLindin = (1:size(dyind,1)).'; 
      dyLindin = dyLindin(nzlocsyL);
      dyLref   = cadaindprint(dyLindin);
      TD2 = ['cada',NDstr,'td2'];
      if xvec
        fprintf(fid,[indent,TD2,' = ',y.deriv(Vcount).name,'(:',dyLref,');']);
      else
        fprintf(fid,[indent,TD2,' = ',y.deriv(Vcount).name,'(',dyLref,');']);
      end
      dyLstr = TD2;
    else
      dyLstr = y.deriv(Vcount).name;
    end
    % Get xL string
    xLref = [jref;kref;iref];
    if toughref; xLref = xLref(tranref); end
    xLref = nonzeros(diag(xLref)*logical(dyL));
    xLref = cadaindprint(xLref(dyLinds));
    if xvec == 1
      xLstr = [x.func.name,'(:,',xLref,')'];
    elseif xvec == 2
      xLstr = [x.func.name,'(',xLref,',:).'''];
    elseif xMrow == 1
      xLstr = [x.func.name,'(',xLref,').'''];
    else
      xLstr = [x.func.name,'(',xLref,')'];
    end
    % Get assignment to dz
    dzyLind  = full(dz(logical(dyL)));
    dzyLasgn = cadaindprint(dzyLind(dyLinds));
    if xvec
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(size(',y.deriv(Vcount).name,',1),%1.0f);\n'],nz);
      end
      dzyLstr  = [TD1,'(:,',dzyLasgn,')'];
    else
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nz);
      end
      dzyLstr  = [TD1,'(',dzyLasgn,')'];
    end
    % Print out the calculation
    if ~dzinit
      fprintf(fid,[indent,dzyLstr,' = ',dyLstr,'.*',xLstr,';\n']);
      dzinit = 1;
    else
      fprintf(fid,[indent,dzyLstr,' = ',dzyLstr,' + ',dyLstr,'.*',xLstr,';\n']);
    end
  end
  %% --------------------- Print out dyR Calcs ------------------------- %%
  if DPFLAG && ~isempty(spdyRind) % dyR
    [~,dyRinds] = sort(nonzeros(dyR));
    % Get dyR string
    if spdyRflag
      dyRindin = (1:size(dyind,1)).'; 
      dyRindin = dyRindin(nzlocsyR);
      dyRref   = cadaindprint(dyRindin);
      TD2 = ['cada',NDstr,'td2'];
      if xvec
        fprintf(fid,[indent,TD2,' = ',y.deriv(Vcount).name,'(:',dyRref,');']);
      else
        fprintf(fid,[indent,TD2,' = ',y.deriv(Vcount).name,'(',dyRref,');']);
      end
      dyRstr = TD2;
    else
      dyRstr = y.deriv(Vcount).name;
    end
    % Get xR string
    xRref = [kref;iref;jref];
    if toughref; xRref = xRref(tranref); end
    xRref = nonzeros(diag(xRref)*logical(dyR));
    xRref = cadaindprint(xRref(dyRinds));
    if xvec == 1
      xRstr = [x.func.name,'(:,',xRref,')'];
    elseif xvec == 2
      xRstr = [x.func.name,'(',xRref,',:).'''];
    elseif xMrow == 1
      xRstr = [x.func.name,'(',xRref,').'''];
    else
      xRstr = [x.func.name,'(',xRref,')'];
    end
    % Get assignment to dz
    dzyRind  = full(dz(logical(dyR)));
    dzyRasgn = cadaindprint(dzyRind(dyRinds));
    if xvec
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(size(',y.deriv(Vcount).name,',1),%1.0f);\n'],nz);
      end
      dzyRstr  = [TD1,'(:,',dzyRasgn,')'];
    else
      if ~dzinit
        fprintf(fid,[indent,TD1,' = zeros(%1.0f,1);\n'],nz);
      end
      dzyRstr  = [TD1,'(',dzyRasgn,')'];
    end
    % Print out the calculation
    if ~dzinit
      fprintf(fid,[indent,dzyRstr,' = -',dyRstr,'.*',xRstr,';\n']);
      dzinit = 1;
    else
      fprintf(fid,[indent,dzyRstr,' = ',dzyRstr,' - ',dyRstr,'.*',xRstr,';\n']);
    end
  end
  if DPFLAG && dzflag
    % Assign to derivstr
    fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
  end
end

if PFLAG
  if dumbflag
    fprintf(fid,[indent,funcstr,' = cross(',x.func.name,',',y.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = cross(',x.func.name,',',y.func.name,',%1.0f);\n'],dim);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id],1) = ADIGATOR.VARINFO.COUNT;
z = class(z,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
