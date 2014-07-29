function yi = ppval(pp,xi)
% CADA overloaded PPVAL
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
if isa(pp,'cada')
  % This came from overloaded interp1
  if ADIGATOR.EMPTYFLAG 
    yi = cadaEmptyEval(pp,xi);
    return
  end
  ADIGATOR.VARINFO.LASTOCC(pp.id,1) = ADIGATOR.VARINFO.COUNT;
  fpp = pp.func.value;
  if ~isfield(fpp,'form')
    fpp.form = 'pp';
  end
  ppknown = pp.func.pp;
elseif isstruct(pp) 
  if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL
    dummybl = length(pp.breaks);
    dummycl = size(pp.coefs,1);
    ADIGATOR.VARINFO.LASTOCC([dummybl.id dummycl.id],1) = ADIGATOR.VARINFO.COUNT;
    if PFLAG
      forOverFlag = isempty(dummybl.func.value);
    else
      forOverFlag = 0;
    end
  end
  if ADIGATOR.EMPTYFLAG
    yi = cadaEmptyEval(pp.breaks,pp.coefs,pp.pieces,pp.order,pp.dim,xi);
    return
  end
  fpp = pp;
  ppknown = 1;
%   if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL
%     OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.breaks.id,:));
%     if length(OverID) > 1; OverID = OverID(1); end
%     if ~PFLAG
%       if ~isempty(OverID)
%         ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
%       end
%       OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.coefs.id,:));
%       if ~isempty(OverID)
%         if length(OverID) > 1; OverID = OverID(1); end
%         ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
%       end
%       OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.pieces.id,:));
%       if ~isempty(OverID)
%         if length(OverID) > 1; OverID = OverID(1); end
%         ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
%       end
%     elseif PFLAG && isempty(OverID)
%       ppknown = 0;
%     end
%   end
  
  
  if ~isfield(fpp,'form')
    fpp.form = 'pp';
  end
  ppnew.id = [];
  ppnew.func = struct('name',[],'size',[1 1],'zerolocs',[],'value',pp,'pp',[]);
  
  fpp.breaks = fpp.breaks.func.value;
  fpp.coefs  = fpp.coefs.func.value;
  fpp.pieces = fpp.pieces.func.value;
  fpp.order  = fpp.order.func.value;
  fpp.dim    = fpp.dim.func.value;
  if (isempty(fpp.breaks) && prod(pp.breaks.func.size)>1) || ...
      (isempty(fpp.coefs) && prod(pp.coefs.func.size)>1) || ...
      (isempty(fpp.pieces) && prod(pp.pieces.func.size)>1)
    % The breaks, coefs, or pieces are unknown..
    ppknown = 0;
    if (isempty(fpp.order) && prod(pp.order.func.size)>1) || ...
        (isempty(fpp.dim) && prod(pp.dim.func.size)>1)
      if ADIGATOR.FORINFO.FLAG && PFLAG
        error(['If using ppval within a loop which you wish to keep rolled, ',...
          'then the pp order and dim must stay consistent across all calls']);
      else
        error('The pp.order and pp.dim must be known numeric inputs')
      end
    end
  end
  if ~ppknown && ~(ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && PFLAG)
    error('all values of piecewise polys must be known numerically')
  end
  ADIGATOR.VARINFO.LASTOCC([pp.breaks.id pp.coefs.id...
    pp.pieces.id pp.order.id pp.dim.id],1) = ADIGATOR.VARINFO.COUNT;
  if PFLAG
    nameloc = strfind(pp.breaks.func.name,'.breaks');
    if ~isempty(nameloc)
      ppnew.func.name = pp.breaks.func.name(1:nameloc(end)-1);
    else
      error('Cannot find the name of the pp - please report error to sourceforge forums');
    end
  end
  pp = ppnew;
else
  error ('invalid input to ppval');
end

if isa(xi,'cada')
  xiMrow = xi.func.size(1);
  xiNcol = xi.func.size(2);
elseif isnumeric(xi)
  [xiMrow, xiNcol] = size(xi);
  xinew.id = [];
  xinew.func = struct('name',[],'size',[xiMrow,xiNcol],'zerolocs',[],'value',xi);
  if PFLAG
    xinew.func.name = cadamatprint(xi);
  end
  xinew.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  xi = xinew;
else
  error('invalid input to ppval')
end

% --------------------------- Build Function ---------------------------- %
if xiMrow > 1 && xiNcol > 1 && fpp.dim > 1
  error('Result of ppval must be vector or matrix')
elseif fpp.dim == 1
  yiMrow = xiMrow;
  yiNcol = xiNcol;
elseif xiMrow == 1 && xiNcol == 1
  yiMrow = fpp.dim;
  yiNcol = 1;
else
  yiMrow = xiMrow*xiNcol;
  yiNcol = fpp.dim;
end
yi.id = ADIGATOR.VARINFO.COUNT;
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
if ~isempty(xi.func.value) && ppknown
  yi.func.value = ppval(fpp,xi.func.value);
end
if isinf(yiMrow); 
  yvec = 1; yiMrow = 1; 
elseif isinf(yiNcol)
  yvec = 2; yiNcol = 1;
else
  yvec = 0;
end
if isinf(xiMrow); xiMrow = 1; elseif isinf(xiNcol); xiNcol = 1; end

if forOverFlag
  ppStr = pp.func.name;
  dblStr = dummybl.func.name;
  dclStr = dummycl.func.name;
  fppStr = ['cada',NDstr,'fpp'];
  fprintf(fid,[indent,fppStr,' = struct(''form'',''pp'',''breaks'',',...
    ppStr,'.breaks(1:',dblStr,'),''coefs'',',...
    ppStr,'.coefs(1:',dclStr,',:),''order'',%1.0d,''pieces'',',...
    ppStr,'.pieces,''dim'',%1.0d);\n'],fpp.order,fpp.dim);
  pp.func.name = fppStr;
end

% --------------------------- Build Derivative -------------------------- %
yi.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if fpp.order > 1 && cadaCheckForDerivs(xi)
  % Build derivative pp
  if DPFLAG && ppknown
    dpp = fpp;
    dpp.coefs(:,end) = [];
    for J = 1:fpp.order-2
      dpp.coefs(:,J) = (fpp.order-J).*dpp.coefs(:,J);
    end
    dpp.order = dpp.order-1;
    TF1 = ['cada',NDstr,'tf1'];
    dppStr = cadamatprint(dpp);
    fprintf(fid,[indent,TF1,' = ppval(',dppStr,',',xi.func.name,');\n']);
  elseif DPFLAG
    dppStr = ['cada',NDstr,'dpp'];
    fprintf(fid,[indent,dppStr,'.form = ''pp'';\n']);
    fprintf(fid,[indent,dppStr,'.breaks = ',pp.func.name,'.breaks;\n']);
    fprintf(fid,[indent,dppStr,'.coefs = [']);
    for J = 1:fpp.order-1
      if fpp.order-J > 1
        fprintf(fid,['%1.0d.*',pp.func.name,'.coefs(:,%1.0d) '],fpp.order-J,J);
      else
        fprintf(fid,[pp.func.name,'.coefs(:,%1.0d) '],J);
      end
    end
    fprintf(fid,'];\n');
    fprintf(fid,[indent,dppStr,'.pieces = ',pp.func.name,'.pieces;\n']);
    fprintf(fid,[indent,dppStr,'.order = %1.0d;\n'],fpp.order-1);
    fprintf(fid,[indent,dppStr,'.dim = %1.0d;\n'],fpp.dim);
    
    TF1 = ['cada',NDstr,'tf1'];
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
  fprintf(fid,[indent,funcstr,' = ppval(',pp.func.name,',',xi.func.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([xi.id yi.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
yi = class(yi,'cada');
return
end