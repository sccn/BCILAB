function y = cadaunarylogical(x,callerstr,dim)
% Called by: logical, not, any, all, isempty
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;

if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end

% ---------------------- Build Function Properties -----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
y.func = struct('name',funcstr,'size',[],'zerolocs',...
  [],'value',[],'logical',[]);
% Function Size
if dim == 1
  y.func.size = [1 x.func.size(2)];
elseif dim == 2
  y.func.size = [x.func.size(1) 1];
elseif dim == -1
  y.func.size = [1 1];
else
  y.func.size = x.func.size;
end

% Function Value/Zero Locations
switch callerstr
  case 'logical'
    if ~isempty(x.func.value)
      y.func.value = logical(x.func.value);
    elseif ~isempty(x.func.zerolocs)
      y.func.zerolocs = x.func.zerolocs;
    end
  case 'not'
    if ~isempty(x.func.value)
      y.func.value = ~logical(x.func.value);
    end
  case 'any'
    if ~isempty(x.func.value)
      y.func.value = any(x.func.value,dim);
    end
  case 'all'
    if ~isempty(x.func.value)
      y.func.value = logical(x.func.value);
    elseif ~isempty(x.func.zerolocs)
      x.func.size(isinf(x.func.size)) = 1;
      xtemp = true(x.func.size);
      xtemp(x.func.zerolocs) = false;
      ytemp = all(xtemp,dim);
      y.func.zerolocs = find(~ytemp(:));
    end
  otherwise
    error(['case not coded: ',callerstr])
end

if PFLAG
  if dim
    fprintf(fid,[indent,funcstr,' = ',...
      callerstr,'(',x.func.name,',%1.0d);\n'],dim);
  else
    fprintf(fid,[indent,funcstr,' = ',...
      callerstr,'(',x.func.name,');\n']);
  end
end

y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = class(y,'cada');
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
return
end