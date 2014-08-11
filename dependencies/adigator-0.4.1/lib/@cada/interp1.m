function yi = interp1(varargin)
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
% -------------------------- Parse Inputs ------------------------------- %
ppflag = 0;
if nargin == 2 
  % yi = interp1(y,xi)
  [y,yMrow,yNcol] = parseinput(varargin{1});
  [xi,xiMrow,xiNcol] = parseinput(varargin{2});
  x.id = []; xMrow = 1; xNcol = yMrow*yNcol;
  x.func = struct('name',[],'size',[xMrow xNcol],...
    'zerolocs',[],'value',1:xNcol); 
  x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  method = 'linear';
elseif nargin == 3
  % yi = interp1(x,y,xi)
  [x,xMrow,xNcol] = parseinput(varargin{1});
  [y,yMrow,yNcol] = parseinput(varargin{2});
  [xi,xiMrow,xiNcol] = parseinput(varargin{3});
  method = 'linear';
elseif nargin == 4 && ischar(varargin{3}) && ...
    ischar(varargin{4}) && strcmp(varargin{4},'pp')
  % pp = interp1(x,y,method,'pp')
  ppflag = 1;
  [x,xMrow,xNcol] = parseinput(varargin{1});
  [y,yMrow,yNcol] = parseinput(varargin{2});
  method = varargin{3};
elseif nargin == 4 && ischar(varargin{3})
  % yi = interp1(y,xi,method,'extrap') or
  % yi = interp1(y,xi,method,extrapval)
  [y,yMrow,yNcol] = parseinput(varargin{1});
  [xi,xiMrow,xiNcol] = parseinput(varargin{2});
  x.id = []; xMrow = 1; xNcol = yMrow*yNcol;
  x.func = struct('name',[],'size',[xMrow xNcol],...
    'zerolocs',[],'value',1:xNcol);
  x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  method    = varargin{3};
  extrapval = varargin{4};
elseif nargin == 4 && ~ischar(varargin{3}) && ischar(varargin{4})
  % yi = interp1(x,y,xi,method)
  [x,xMrow,xNcol] = parseinput(varargin{1});
  [y,yMrow,yNcol] = parseinput(varargin{2});
  [xi,xiMrow,xiNcol] = parseinput(varargin{3});
  method = varargin{4};
elseif nargin == 5
  % yi = interp1(x,y,xi,method,'extrap') or
  % yi = interp1(x,y,xi,method,extrapval)
  [x,xMrow,xNcol] = parseinput(varargin{1});
  [y,yMrow,yNcol] = parseinput(varargin{2});
  [xi,xiMrow,xiNcol] = parseinput(varargin{3});
  method = varargin{4};
  extrapval = varargin{5};
  if ~ischar(extrapval)
    warning('CADA overloaded interp1 does not allow for extrapvals to be given'); %#ok<WNTAG>
  end
else
  error('Invalid input to interp1')
end
if ADIGATOR.EMPTYFLAG
  if ppflag
    yi = cadaEmptyEval(x,y);
  else
    yi = cadaEmptyEval(x,y,xi);
  end
  return
end
if ~exist('extrapval','var')
  switch method(1)
    case {'s','p','c'}
      extrapval = 'extrap';
    otherwise
      extrapval = NaN;
  end
end
% Check x/y
if isinf(xMrow) || isinf(xNcol) || isinf(yMrow) || isinf(yNcol)
  error('X and Y must be non-vectorized')
elseif isempty(x.func.value) || isempty(y.func.value)
  if ADIGATOR.FORINFO.FLAG && PFLAG
    errstr = ['Cannot use interp1 within a loop if the values of ', ...
      'X and Y change on each iteration - please generate the piecewise ',...
      'polynomials using interp1, spline, etc. PRIOR to differentiation and ',...
      'input them as known numeric inputs, ppval will then allow for you to ',...
      'loop on the different generate polys- see manual for more information'];
    error(errstr)
  else
    error('X and Y inputs to interp1 must be known numerically')
  end
elseif xMrow > 1 && xNcol > 1
  error('X must be a vector')
elseif cadaCheckForDerivs(x) || cadaCheckForDerivs(y)
  error('Can only take derivative of interp1 wrt input XI, not X or Y')
end

yi.id = ADIGATOR.VARINFO.COUNT;
% -------------------- Piecewise Polynomial Output ---------------------- %
if ppflag
  % We are going to use some duct tape.. Going to make the output a cada
  % object
  funcstr = cadafuncname();
  fpp = interp1(x.func.value,y.func.value,method,'pp');
  yi.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',fpp,'pp',ppknown);
  if PFLAG
    fprintf(fid,[indent,funcstr,' = ',cadamatprint(fpp),';\n']);
  end
  yi.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  ADIGATOR.VARINFO.LASTOCC([x.id y.id yi.id],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  yi = class(yi,'cada');
  return
elseif PFLAG
  fpp = interp1(x.func.value,y.func.value,method,'pp');
end

% ------------------------- Build Function ------------------------------ %
if xiMrow > 1 && xiNcol > 1 && yMrow > 1 && yNcol > 1
  error('Overloaded interp1 must return a vector or matrix')
elseif yMrow > 1 && yNcol > 1
  if yMrow ~= xMrow*xNcol
    error('length(X) and size(Y,1) must be the same')
  end
  yiMrow = xiMrow*xiNcol;
  yiNcol = yNcol;
  if isinf(xiMrow) || isinf(xiNcol)
    error(['if XI is vectorized, output of interp1 must be N by n, ',...
    'or n by N, where N is length of vectorized dim'])
  end
elseif xMrow*xNcol == yMrow*yNcol
  yiMrow = xiMrow;
  yiNcol = xiNcol;
else
  error('X and Y must be of same length')
end

[funcstr,DPFLAG] = cadafuncname();
yi.func = struct('name',funcstr,'size',[yiMrow yiNcol],'zerolocs',...
  [],'value',[]);
if isfield(xi.func,'logicref')
  if fpp.dim == 1
    yi.func.logicref = xi.func.logicref;
  elseif xiMrow == 1 && ~xi.func.logicref(1)
    yi.func.logicref = [xi.func.logicref(1), 0];
  elseif xiNcol == 1 && ~xi.func.logicref(2)
    yi.func.logicref = xi.func.logicref;
  else
    error('Cannot track the logical reference in this instance');
  end
end
if ~isempty(xi.func.value)
  yi.func.value = interp1(x.func.value,y.func.value,xi.func.value,method,extrapval);
end
if isinf(yiMrow); 
  yvec = 1; yiMrow = 1; 
elseif isinf(yiNcol)
  yvec = 2; yiNcol = 1;
else
  yvec = 0;
end
if isinf(xiMrow); xiMrow = 1; elseif isinf(xiNcol); xiNcol = 1; end

if PFLAG
  ppStr = cadamatprint(fpp);
end

% --------------------------- Build Derivative -------------------------- %
yi.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if fpp.order > 1 && cadaCheckForDerivs(xi)
  % Build derivative pp
  if DPFLAG
    dpp = fpp;
    dpp.coefs(:,end) = [];
    for J = 1:fpp.order-2
      dpp.coefs(:,J) = (fpp.order-J).*dpp.coefs(:,J);
    end
    dpp.order = dpp.order-1;
    TF1 = ['cada',NDstr,'tf1'];
    dppStr = cadamatprint(dpp);
    fprintf(fid,[indent,TF1,' = ppval(',dppStr,',',xi.func.name,');\n']);
  end
  for Vcount = 1:NUMvod
    if ~isempty(xi.deriv(Vcount).nzlocs)
      derivstr = cadadername(funcstr,Vcount);
      yi.deriv(Vcount).name    = derivstr;
      if xiMrow == yiMrow && xiNcol == yiNcol
        % yi is same size as xi
        yi.deriv(Vcount).nzlocs = xi.deriv(Vcount).nzlocs;
        DYiStr = xi.deriv(Vcount).name;
        yirows = yi.deriv(Vcount).nzlocs(:,1);
      else
        % need to repmat xi derivs first
        xirows = xi.deriv(Vcount).nzlocs(:,1);
        xicols = xi.deriv(Vcount).nzlocs(:,2);
        nzxi   = length(xirows);
        nv     = ADIGATOR.VAROFDIFF(Vcount).usize;
        dxi = sparse(xirows,xicols,1:nzxi,xiMrow*xiNcol,nv);
        dyi = repmat(dxi,[yiNcol 1]);
        [yirows, yicols, xilocs] = find(dyi);
        if size(yirows,2) > 1
          yirows = yirows.'; yicols = yicols.';
        end
        yi.deriv(Vcount).nzlocs = [yirows yicols];
        if DPFLAG
          DYiStr  = ['cada',NDstr,'td1'];
          xindstr = cadaindprint(xilocs);
          if yvec
            fprintf(fid,[indent,DYiStr,' = ',x.deriv(Vcount).name,'(:,',xindstr,');\n']);
          else
            fprintf(fid,[indent,DYiStr,' = ',x.deriv(Vcount).name,'(',xindstr,');\n']);
          end
        end
      end
      if DPFLAG
        TF2 = ['cada',NDstr,'tf2'];
        findstr = cadaindprint(yirows);
        if yvec == 1
          fprintf(fid,[indent,derivstr,' = ',TF1,'(:,',findstr,').*',DYiStr,';\n']);
        elseif yvec
          fprintf(fid,[indent,derivstr,' = ',TF1,'(',findstr,',:).''.*',DYiStr,';\n']);
        else
          fprintf(fid,[indent,TF2,' = ',TF1,'(',findstr,');\n']);
          fprintf(fid,[indent,derivstr,' = ',TF2,'(:).*',DYiStr,';\n']);
        end
      end
    end
  end
end

% Print Function
if PFLAG
  if xiMrow == 1 && xiNcol == 1 && yMrow > 1 && yNcol > 1
    fprintf(fid,[indent,funcstr,' = ppval(',ppStr,',',xi.func.name,').'';\n']);
  else
    fprintf(fid,[indent,funcstr,' = ppval(',ppStr,',',xi.func.name,');\n']);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id yi.id xi.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
yi = class(yi,'cada');
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