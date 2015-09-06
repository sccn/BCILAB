function y = nonzeros(x)
% CADA overloaded version of function NONZEROS.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);

if ADIGATOR.FORINFO.FLAG
  IncreaseForNzCount();
  if ADIGATOR.RUNFLAG == 2
    y = ForNonzeros(x);
    return
  end
end
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end
% --------------------Build Function Properties---------------------------%
xMrow = x.func.size(1); xNcol = x.func.size(2);
if isinf(xMrow) || isinf(xNcol)
  error('Cannot use non-zeros on vectorized object')
end
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[],'zerolocs',[],'value',[]);

y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if ~isempty(x.func.value)
  % ---x has numeric values--- %
  y.func.value  = nonzeros(x.func.value);
  y.func.size   = size(y.func.value);
  % - no derivatives
  % Function Printing
  if PFLAG
    fprintf(fid,[indent,funcstr,' = nonzeros(',x.func.name,');\n']);
  end
  if ADIGATOR.FORINFO.FLAG
    AssgnForNzInds(find(x.func.value(:)),1,x.func.size.');
  end
elseif ~isempty(x.func.zerolocs)
  % ---x has zero entries--- %
  isubs = true(xMrow, xNcol);   isubs(x.func.zerolocs) = false;
  isubs = find(isubs(:));
  y.func.size = [length(isubs),1];
  % ------------------Build Derivative Properties-------------------------%
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).name)
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
      if DPFLAG == 1
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
    end
  end
  if ADIGATOR.FORINFO.FLAG
    AssgnForNzInds(isubs,3,x.func.size.');
  end
  if PFLAG == 1
    TF1 = ['cada',NDstr,'tf1'];
    TFind1 = cadaindprint(isubs);
    fprintf(fid,[TF1,' = zeros(%1.0f,1);\n'],length(isubs));
    fprintf(fid,[TF1,'(:) = ',x.func.name,'(',TFind1,');\n']);
    fprintf(fid,[funcstr,' = ',TF1,';\n']);
  end
else
  % ---x is full - y is just x reshaped--- %
  y.func.size = [xMrow*xNcol 1];
  % ------------------Build Derivative Properties-------------------------%
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).name)
      derivstr = cadadername(funcstr,Vcount);
      y.deriv(Vcount).name = derivstr;
      y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
      if DPFLAG == 1
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
    end
  end
  if PFLAG == 1
    fprintf(fid,[indent,funcstr,' = ',x.func.name,'(:);\n']);
  end
  if ADIGATOR.FORINFO.FLAG
    AssgnForNzInds(1:xMrow*xNcol,2,x.func.size.');
  end
end

ADIGATOR.FUNCTIONLOCATION([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
return
end

function IncreaseForNzCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.NONZEROS =...
  ADIGATORFORDATA(INNERLOC).COUNT.NONZEROS + 1;
return
end

function AssgnForNzInds(inds,flag,xsize)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
NZCOUNT   = ADIGATORFORDATA(INNERLOC).COUNT.NONZEROS;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
if isempty(ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS)
  ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).INDICES  = inds(:);
  ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS    = [flag length(inds)];
else
  if ITERCOUNT == 1
    ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).INDICES  = inds(:);
  elseif ~isempty(inds)
    ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).INDICES(1:length(inds),ITERCOUNT) = inds;
  end
  pflag = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(1);
  yRow  = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(2);
  
  ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(1) = max(flag,pflag);
  ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(2) = max(yRow,length(inds));
end

ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).SIZES = ...
  [ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).SIZES,xsize];
return
end

function y = ForNonzeros(x)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
NZCOUNT   = ADIGATORFORDATA(INNERLOC).COUNT.NONZEROS;
flag      = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(1);
yMrow     = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).FLAGS(2);
NUMvod    = ADIGATOR.NVAROFDIFF;
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
DNstr     = sprintf('%1.0d',ADIGATOR.DERNUMBER);


y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func  = struct('name',funcstr,'size',[yMrow 1],'zerolocs',[],'value',[]);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
% ---Print Out Derivatives--- %
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    y.deriv(Vcount).name = derivstr;
    y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
    if DPFLAG == 1
      fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
    end
  end
end

% ---Print Out Functions---%
IndName     = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).INDICES{1};
tempfuncstr = ['cada',DNstr,'tempf1'];
if ~isempty(IndName)
  % Need to use references which have already been printed.
  CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;
  IndDep    = ADIGATORFORDATA(INNERLOC).NONZEROS(NZCOUNT).INDICES{3}(1);
  if IndDep
    fprintf(fid,[indent,tempfuncstr,' = ',x.func.name,'(nonzeros(',IndName,'(:,',CountName,')));\n']);
  else
    fprintf(fid,[indent,tempfuncstr,' = ',x.func.name,'(nonzeros(',IndName,'));\n']);
  end
elseif flag == 1
  fprintf(fid,[indent,tempfuncstr,' = nonzeros(',x.func.name,');\n']);
elseif flag == 2
  fprintf(fid,[indent,tempfuncstr,' = ',x.func.name,'(:);\n']);
end
fprintf(fid,[indent,y.func.name,' = zeros(%1.0d,1);\n'],y.func.size(1));
fprintf(fid,[indent,y.func.name,'(1:length(',tempfuncstr,')) = ',tempfuncstr,';\n']);

ADIGATOR.FUNCTIONLOCATION([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
return
end