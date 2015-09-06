function y = full(x)
% CADA overloaded FULL function.
%
% Written by Matthew J. Weinstein
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
% --------------------Build Function Properties-----------------------%
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',x.func.size,'zerolocs',...
  x.func.zerolocs,'value',x.func.value);

% ------------------Build Derivative Properties-----------------------%
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
        derivstr = cadadername(funcstr,Vcount);
        y.deriv(Vcount).name    = derivstr;
        y.deriv(Vcount).nzlocs  = x.deriv(Vcount).nzlocs;
        % -----------------------Derivative Printing----------------------%
        if DPFLAG == 1
          fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
        end
    end
end

% --------------------------Function Printing ----------------------------%
if PFLAG == 1
    fprintf(fid,[indent,funcstr,' = full(',x.func.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = class(y,'cada');
