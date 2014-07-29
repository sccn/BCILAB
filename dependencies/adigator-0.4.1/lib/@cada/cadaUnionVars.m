function z = cadaUnionVars(x,y)
% Union two CADA variables to unite sizes and sparsity patterns
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ~isa(x,'cada')
  z = y;
  return
elseif ~isa(y,'cada')
  z = x;
  return
end
z.id = max([x.id y.id]);
% Function Sizing
xMrow = x.func.size(1); xNcol = x.func.size(2);
yMrow = y.func.size(1); yNcol = y.func.size(2);
% Vectorized Stuff
vecerr = 0;
if (isinf(xMrow) && isinf(yMrow)) || (isinf(xMrow) && ~yMrow) ...
    || (~xMrow && isinf(yMrow))
  if isinf(xNcol) || isinf(yNcol); vecerr=1; end
  zvec = 1;
elseif (isinf(xNcol) && isinf(yNcol)) || (isinf(xNcol) && ~yNcol) ...
    || (~xNcol && isinf(yNcol))
  if isinf(xMrow) || isinf(yMrow); vecerr=1; end
  zvec = 2;
else
  if isinf(xNcol) || isinf(yNcol) || isinf(xMrow) || isinf(yMrow)
    vecerr = 1;
  end
  zvec = 0;
end
if vecerr
  error(['If a variable may take on different sizes depending upon ',...
    'loop iterations or different branches of a conditional block ',...
    'then it must not switch between vectorized to non-vectorized'])
end
if xMrow > yMrow; zMrow = xMrow; else zMrow = yMrow; end
if xNcol > yNcol; zNcol = xNcol; else zNcol = yNcol; end

z.func = struct('name',[],'size',[zMrow zNcol],'zerolocs',[],'value',[]);
if zvec == 1; 
  xMrow = 1; yMrow = 1; zMrow = 1; 
elseif zvec == 2; 
  xNcol = 1; yNcol = 1; zNcol = 1; 
end


% Function Numerics/Sparsity
if isfield(x.func,'pp') && isfield(y.func,'pp')
  % Overloaded piecewise polys
  if x.func.pp && y.func.pp && isequal(x.func.value,y.func.value)
    z.func.pp = 1;
    z.func.value = x.func.value;
  elseif strcmp(x.func.value.form,'pp') && strcmp(y.func.value.form,'pp')
    if x.func.value.order ~= y.func.value.order
      error('Attempting to union two peicewise polys where the order is different - cannot do this')
    elseif x.func.value.dim ~= y.func.value.dim
      error('Attempting to union two peicewise polys where the dim is different - cannot do this')
    end
    z.func.value  = x.func.value;
    z.func.breaks = [];
    z.func.coefs  = [];
    z.func.pieces = [];
    z.func.pp = 0;
  elseif strcmp(x.func.value.form,'adigatorpp2') && strcmp(y.func.value.form,'adigatorpp2')
    if x.func.value.xorder ~= y.func.value.xorder
      error('Attempting to union two 2D peicewise polys where the xorder is different - cannot do this')
    elseif x.func.value.yorder ~= y.func.value.yorder
      error('Attempting to union two 2D peicewise polys where the yorder is different - cannot do this')
    end
    z.func.value  = x.func.value;
    z.func.xbreaks = []; z.func.ybreaks = [];
    z.func.coefs  = [];
    z.func.pp = 0;
  else
    error('Cannot determine types of polynomials attempting to union -  please report bug')
  end
elseif ~isempty(x.func.value) && ~isempty(y.func.value)
  if isequal(x.func.value,y.func.value)
    z.func.value = x.func.value;
  else
    xtemp = logical(x.func.value); ytemp = logical(y.func.value);
  end
elseif ~isempty(x.func.value)
  xtemp = logical(x.func.value);
  ytemp = true(yMrow,yNcol);
  if ~isempty(y.func.zerolocs); ytemp(y.func.zerolocs) = false; end
elseif ~isempty(y.func.value)
  ytemp = logical(y.func.value);
  xtemp = true(xMrow,xNcol);
  if ~isempty(x.func.zerolocs); xtemp(x.func.zerolocs) = false; end
else
  xtemp = true(xMrow,xNcol); ytemp = true(yMrow,yNcol);
  if ~isempty(x.func.zerolocs); xtemp(x.func.zerolocs) = false; end
  if ~isempty(y.func.zerolocs); ytemp(y.func.zerolocs) = false; end
end
if isempty(z.func.value)
  if zMrow > xMrow || zNcol > xNcol
    xtemp(zMrow,zNcol) = false;
  end
  if zMrow > yMrow || zNcol > yNcol
    ytemp(zMrow,zNcol) = false;
  end
  ztemp = xtemp+ytemp;
  z.func.zerolocs = find(~ztemp(:));
  if length(z.func.zerolocs) == zMrow*zNcol
    z.func.zerolocs = []; z.func.value = zeros(zMrow,zNcol);
  end
end

if (isfield(x.func,'logicref') && ~isfield(y.func,'logicref')) || ...
    (~isfield(x.func,'logicref') && isfield(y.func,'logicref'))
  error('Attempting to union one variable which results from a logical reference with one which does not');
elseif isfield(x.func,'logicref')
  if ~isequal(x.func.logicref,y.func.logicref)
    error('Cannot union two variables if they result from different logical references');
  else
    z.func.logicref = x.func.logicref;
  end
end

% Derivative Sparsity
NUMvod = ADIGATOR.NVAROFDIFF;
z.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs) && ~isempty(y.deriv(Vcount).nzlocs)
    z.deriv(Vcount).name = y.deriv(Vcount).name;
    nv    = ADIGATOR.VAROFDIFF(Vcount).usize;
    xrows = x.deriv(Vcount).nzlocs(:,1);
    xcols = x.deriv(Vcount).nzlocs(:,2);
    yrows = y.deriv(Vcount).nzlocs(:,1);
    ycols = y.deriv(Vcount).nzlocs(:,2);
    if xMrow~=yMrow &&zNcol > 1
      zrows = [xrows + (zMrow - xMrow)*floor((xrows-1)/xMrow);...
        yrows + (zMrow - yMrow)*floor((yrows-1)/yMrow)];
    else
      zrows = [xrows;yrows];
    end
    zcols = [xcols;ycols];
    ztemp = sparse(zrows,zcols,ones(size(zrows)),zMrow*zNcol,nv);
    [zrows,zcols] = find(ztemp);
    if size(zrows,2) > 1
      zrows = zrows.'; zcols = zcols.';
    end
    z.deriv(Vcount).nzlocs = [zrows,zcols];
  elseif ~isempty(x.deriv(Vcount).nzlocs)
    z.deriv(Vcount) = x.deriv(Vcount);
    if xMrow~=yMrow && zNcol > 1
      xrows = x.deriv(Vcount).nzlocs(:,1); xcols = x.deriv(Vcount).nzlocs(:,2);
      zrows = xrows + (zMrow - xMrow)*floor((xrows-1)/xMrow);
      z.deriv(Vcount).nzlocs = [zrows,xcols];
    end
  elseif ~isempty(y.deriv(Vcount).nzlocs)
    z.deriv(Vcount) = y.deriv(Vcount);
    if xMrow~=yMrow && zNcol > 1
      yrows = y.deriv(Vcount).nzlocs(:,1); 
      ycols = y.deriv(Vcount).nzlocs(:,2);
      zrows = yrows + (zMrow - yMrow)*floor((yrows-1)/yMrow);
      z.deriv(Vcount).nzlocs = [zrows,ycols];
    end
  end
end
z = class(z,'cada');
return
end