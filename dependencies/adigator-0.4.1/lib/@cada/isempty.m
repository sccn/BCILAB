function y = isempty(x)
% CADA overloaded ISEMPTY function.
global ADIGATOR
Dbstuff = dbstack; CallingFile = Dbstuff(2).file;
if length(CallingFile) > 12 && ~isempty(strfind(CallingFile,'adigatortempfunc'))
  if ADIGATOR.EMPTYFLAG
    y = cadaEmptyEval(x);
    return
  end
  NUMvod  = ADIGATOR.NVAROFDIFF;
  fid     = ADIGATOR.PRINT.FID;
  PFLAG   = ADIGATOR.PRINT.FLAG;
  indent  = ADIGATOR.PRINT.INDENT;
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  % --------------------Build Function Properties-----------------------%
  y.id = ADIGATOR.VARINFO.COUNT;
  [funcstr,~] = cadafuncname();
  y.func = struct('name',funcstr,'size',[1 1],'zerolocs',...
    [],'value',[],'logical',[]);
  
  if xMrow*xNcol > 0
    y.func.value = false;
  else
    y.func.value = true;
  end
  
  y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  
  % --------------------------Function Printing ----------------------------%
  if PFLAG == 1
    fprintf(fid,[indent,funcstr,' = ~',x.func.name,';\n']);
  end
  
  ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  y = class(y,'cada');
  
else
  y = false;
end