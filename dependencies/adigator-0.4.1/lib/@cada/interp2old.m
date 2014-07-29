function zi = interp2(x,y,z,xi,yi,varargin)
% CADA overloaded version of interp2, note: only the following two input
% structures are allowed:
%       zi = interp2(x,y,z,xi,yi,method);
%                 or
%       zi = interp2(x,y,z,xi,yi);
% Where in the second instance, the method will be treated as the default
% linear. 
%
% NOTE: This will not accept method 'nearest'
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
% -------------------------- Parse Inputs ------------------------------- %
[x, ~, ~] = parseinput(x);
[y, ~, ~] = parseinput(y);
[z, ~, ~] = parseinput(z);
[xi, xiMrow, xiNcol] = parseinput(xi);
[yi, yiMrow, yiNcol] = parseinput(yi);
if ADIGATOR.EMPTYFLAG
  zi = cadaEmptyEval(x,y,z,xi,yi);
  return
end
if nargin == 6
  method = varargin{1};
  if ~strcmp(method,'linear') && ~strcmp(method,'cubic') && ~strcmp(method,'spline')
    error('overloaded interp2 only accepts methods of ''linear'', ''cubic'', or ''spline''');
  end
else
  method = 'linear';
end
% Check X,Y,Z
if isempty(x.func.value) || isempty(y.func.value) || isempty(z.func.value)
  error('X, Y, and Z inputs to interp2 must be known numerically')
end
[msg,X,Y,Z,~] = xyzchk(x.func.value,y.func.value,z.func.value);
error(msg);
% Check XI, YI
if xi.func.size == yi.func.size
  % X, Y, Z same size
  ziMrow = xi.func.size(1); ziNcol = xi.func.size(2);
  repflag = 0;
elseif (xiMrow == 1 && yiNcol == 1) || (xiNcol == 1 && yiMrow == 1)
  % Z takes on row dimension of length of Y and col dimension of length of
  % X
  ziMrow = yiMrow*yiNcol;   ziNcol = xiMrow*xiNcol;
  if isinf(ziMrow) && isinf(ziNcol)
    error('Result of vectorized interp2 may only have 1 vectorized dimension')
  elseif isinf(ziMrow)
    yrep = ones(1,ziNcol);
    if cadaCheckForDerivs(xi)
      error('Cannot perform interp2 when YI vectorized, XI non-vectorized and has derivative information')
    end
  elseif isinf(ziNcol)
    xrep = ones(ziMrow,1);
    if cadaCheckForDerivs(yi)
      error('Cannot perform interp2 when XI vectorized, YI non-vectorized and has derivative information')
    end
  else
    xrep = repmat(1:ziNcol,[ziMrow 1]);     xrep = xrep(:);
    yrep = repmat((1:ziMrow).',[1 ziNcol]); yrep = yrep(:);
  end
  repflag = 1; 
else
  error('XI and YI must be the same size or vectors of different orientations.')
end

% Build ZI func
zi.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
zi.func = struct('name',funcstr,'size',[ziMrow ziNcol],'zerolocs',...
  [],'value',[]);
if ~isempty(xi.func.value) && ~isempty(yi.func.value)
  zi.func.value = interp2(X,Y,Z,xi.func.value,yi.func.value,method);  
end
if isinf(ziMrow); 
  ziMrow = 1; zvec = 1; 
elseif isinf(ziNcol); 
  ziNcol = 1; zvec = 2;
else
  zvec = 0;
end
% Check Logical Reference
if repflag
  if isfield(xi.func,'logicref') || isfield(yi.func,'logicref')
    zi.func.logicref = zeros(1,2);
    if isfield(xi.func,'logicref')
      zi.func.logicref(2) = nonzeros(xi.func.logicref);
    end
    if isfield(yi.func,'logicref')
      zi.func.logicref(1) = nonzeros(yi.func.logicref);
    end
  end
elseif isfield(xi.func,'logicref') && isfield(yi.func,'logicref') && ...
    isequal(xi.func.logicref,yi.func.logicref)
  zi.func.logicref = xi.func.logicref;
elseif isfield(xi.func,'logicref') || isfield(yi.func,'logicref')
  error('Invalid binary operation on result of an unknown logical reference')
end

% Build pp (if printing)
if PFLAG
  pp = adigatorGenInterp2pp(X,Y,Z,method);
end


% Build ZI deriv
ppxflag = 0; ppyflag = 0;
zi.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  nv = ADIGATOR.VAROFDIFF(Vcount).usize;
  % Repmat derivs if needed
  if repflag && ~isempty(xi.deriv(Vcount).nzlocs)
    xrows = xi.deriv(Vcount).nzlocs(:,1);
    xcols = xi.deriv(Vcount).nzlocs(:,2);
    dx = sparse(xrows,xcols,1:length(xrows),ziNcol,nv);
    dx = dx(xrep,:);
    [xrows, xcols,xrepind] = find(dx);
    if size(xrows,2) > 1; xrows = xrows.'; xcols = xcols.'; end
    xi.deriv(Vcount).nzlocs = [xrows, xcols];
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      TDind1 = cadaindprint(xrepind);
      if zvec
        fprintf(fid,[indent,TD1,' = ',xi.deriv(Vcount).name,'(:,',TDind1,');\n']);
      else
        fprintf(fid,[indent,TD1,' = ',xi.deriv(Vcount).name,'(',TDind1,');\n']);
      end
      xi.deriv(Vcount).name = TD1;
    end
  end
  if repflag && ~isempty(yi.deriv(Vcount).nzlocs)
    yrows = yi.deriv(Vcount).nzlocs(:,1);
    ycols = yi.deriv(Vcount).nzlocs(:,2);
    dy = sparse(yrows,ycols,1:length(yrows),ziMrow,nv);
    dy = dy(yrep,:);
    [yrows, ycols,yrepind] = find(dy);
    if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
    yi.deriv(Vcount).nzlocs = [yrows, ycols];
    if DPFLAG
      TD2 = ['cada',NDstr,'td2'];
      TDind1 = cadaindprint(yrepind);
      if zvec
        fprintf(fid,[indent,TD2,' = ',yi.deriv(Vcount).name,'(:,',TDind1,');\n']);
      else
        fprintf(fid,[indent,TD2,' = ',yi.deriv(Vcount).name,'(',TDind1,');\n']);
      end
      yi.deriv(Vcount).name = TD2;
    end
  end
  % Print out dZdX and dZdY if needed
  if DPFLAG && ~ppxflag && ~isempty(xi.deriv(Vcount).nzlocs)
    % Get pp which is for dzi/dxi
    ppx     = getppderiv(pp,'x');
    ppxStr  = cadamatprint(ppx);
    dZdXstr = ['cada',NDstr,'tf1'];
    fprintf(fid,[indent,dZdXstr,' = adigatorEvalInterp2pp(',ppxStr,',',xi.func.name,',',yi.func.name,');\n']);
  end
  if DPFLAG && ~ppyflag && ~isempty(yi.deriv(Vcount).nzlocs)
    ppy     = getppderiv(pp,'y');
    ppyStr  = cadamatprint(ppy);
    dZdYstr = ['cada',NDstr,'tf2'];
    fprintf(fid,[indent,dZdYstr,' = adigatorEvalInterp2pp(',ppyStr,',',xi.func.name,',',yi.func.name,');\n']);
  end
  if ~isempty(xi.deriv(Vcount).nzlocs) && ~isempty(yi.deriv(Vcount).nzlocs)
    % XI and YI have derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    [zi.deriv(Vcount).nzlocs,dzxind,dzyind,xindflag,yindflag] = ...
      cadaunion(xi.deriv(Vcount).nzlocs,yi.deriv(Vcount).nzlocs,...
      ziMrow*ziNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
    xrows = xi.deriv(Vcount).nzlocs(:,1);
    yrows = yi.deriv(Vcount).nzlocs(:,1);
    nzz   = size(zi.deriv(Vcount).nzlocs,1);
    if DPFLAG
      TD3 = ['cada',NDstr,'td3'];
      % Print dZdX
      TF3 = ['cada',NDstr,'tf3'];
      TFind = cadaindprint(xrows);
      if zvec == 1
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(:,',TFind,');\n']);
      elseif zvec == 2
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,',:).'';\n']);
      else
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,');\n']);
      end
      if xindflag
        DZXstr = TD3;
      elseif zvec
        fprintf(fid,[indent,TD3,' = zeros(size(',xi.deriv(Vcount).name,',1),%1.0f);\n'],nzz);
        TDind = cadaindprint(dzxind);
        DZXstr = [TD3,'(:,',TDind,')'];
      else
        fprintf(fid,[indent,TD3,' = zeros(%1.0f,1);\n'],nzz);
        TDind = cadaindprint(dzxind);
        DZXstr = [TD3,'(',TDind,')'];
      end
      if zvec
        fprintf(fid,[indent,DZXstr,' = ',TF3,'.*',xi.deriv(Vcount).name,';\n']);
      else
        fprintf(fid,[indent,DZXstr,' = ',TF3,'(:).*',xi.deriv(Vcount).name,';\n']);
      end
      
      % Print dZdY
      TFind = cadaindprint(yrows);
      if zvec == 1
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(:,',TFind,');\n']);
      elseif zvec == 2
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,',:).'';\n']);
      else
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,');\n']);
      end
      if yindflag
        DZYstr = TD3;
      elseif zvec
        TDind = cadaindprint(dzyind);
        DZYstr = [TD3,'(:,',TDind,')'];
      else
        TDind = cadaindprint(dzyind);
        DZYstr = [TD3,'(',TDind,')'];
      end
      if zvec
        fprintf(fid,[indent,DZYstr,' = ',DZYstr,' + ',TF3,'.*',yi.deriv(Vcount).name,';\n']);
      else
        fprintf(fid,[indent,DZYstr,' = ',DZYstr,' + ',TF3,'(:).*',yi.deriv(Vcount).name,';\n']);
      end
      fprintf(fid,[indent,derivstr,' = ',TD3,';\n']);
    end
  elseif ~isempty(xi.deriv(Vcount).nzlocs)
    % XI has derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    zi.deriv(Vcount).nzlocs = xi.deriv(Vcount).nzlocs;
    if DPFLAG
      % Print dZdX
      xrows = xi.deriv(Vcount).nzlocs(:,1);
      TFind = cadaindprint(xrows);
      if zvec == 1
        fprintf(fid,[indent,derivstr,' = ',dZdXstr,'(:,',TFind,').*',xi.deriv(Vcount).name,';\n']);
      elseif zvec == 2
        fprintf(fid,[indent,derivstr,' = ',dZdXstr,'(',TFind,',:).''.*',xi.deriv(Vcount).name,';\n']);
      else
        TF3 = ['cada',NDstr,'tf3'];
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,');\n']);
        fprintf(fid,[indent,derivstr,' = ',TF3,'(:).*',xi.deriv(Vcount).name,';\n']);
      end
    end
  elseif ~isempty(yi.deriv(Vcount).nzlocs)
    % YI has derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    zi.deriv(Vcount).nzlocs = yi.deriv(Vcount).nzlocs;
    if DPFLAG
      % Print dZdY
      yrows = yi.deriv(Vcount).nzlocs(:,1);
      TFind = cadaindprint(yrows);
      if zvec == 1
        fprintf(fid,[indent,derivstr,' = ',dZdYstr,'(:,',TFind,').*',yi.deriv(Vcount).name,';\n']);
      elseif zvec == 2
        fprintf(fid,[indent,derivstr,' = ',dZdYstr,'(',TFind,',:).''.*',yi.deriv(Vcount).name,';\n']);
      else
        TF3 = ['cada',NDstr,'tf3'];
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,');\n']);
        fprintf(fid,[indent,derivstr,' = ',TF3,'(:).*',yi.deriv(Vcount).name,';\n']);
      end
    end
  end
end

% Print Function
if PFLAG
  ppStr = cadamatprint(pp);
  fprintf(fid,[indent,funcstr,' = adigatorEvalInterp2pp(',ppStr,',',...
    xi.func.name,',',yi.func.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id z.id xi.id yi.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
zi = class(zi,'cada');
return
end

function [y, yMrow, yNcol] = parseinput(x)
if isa(x,'cada')
  y = x; yMrow = x.func.size(1); yNcol = x.func.size(2);
elseif isnumeric(x)
  [yMrow, yNcol] = size(x);
  y.id = [];
  y.func = struct('name',[],'size',[yMrow, yNcol],'zerolocs',[],'value',x);
else
  error('Invalid input to interp1')
end
end

function dpp = getppderiv(pp,dim)
switch dim
  case 'y'
    yorder    = pp.yorder;
    dyorder   = yorder-1;
    coefs     = pp.coefs;
    % Remove Constant Elements
    coefs(:,:,yorder,:) = [];
    % Multiply
    for I = 1:dyorder
      coefs(:,:,I,:) = coefs(:,:,I,:).*(yorder-I);
    end
    % Assign dpp
    dpp = pp;
    dpp.coefs  = coefs;
    dpp.yorder = dyorder;
  case 'x'
    xorder    = pp.xorder;
    dxorder   = xorder-1;
    coefs     = pp.coefs;
    % Remove Constant Elements
    coefs(:,:,:,xorder) = [];
    % Multiply
    for J = 1:dxorder
      coefs(:,:,:,J) = coefs(:,:,:,J).*(xorder-J);
    end
    % Assign dpp
    dpp = pp;
    dpp.coefs  = coefs;
    dpp.xorder = dxorder;
end
end