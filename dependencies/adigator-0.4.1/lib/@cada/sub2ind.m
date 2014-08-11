function K = sub2ind(Asize,I,J)
% Overloaded sub2ind - only works on known numeric objects and with 3
% inputs
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.EMPTYFLAG
  K = cadaEmptyEval(Asize,I,J);
  return
end
PFLAG = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
fid   = ADIGATOR.PRINT.FID;
NUMvod = ADIGATOR.NVAROFDIFF;
[Asize,Astr] = parseInput(Asize);
[I,Istr] = parseInput(I);
[J,Jstr] = parseInput(J);
Knum = sub2ind(Asize,I,J);

K.id = ADIGATOR.VARINFO.COUNT;
[funcstr,~] = cadafuncname();
K.func = struct('name',funcstr,'size',size(Knum),'zerolocs',[],'value',Knum);
K.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if PFLAG
  fprintf(fid,[indent,funcstr,' = sub2ind(',Astr,',',Istr,',',Jstr,');\n']);
end
ADIGATOR.VARINFO.LASTOCC(K.id,1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
K = class(K,'cada');
end
function [num,str] = parseInput(x)
global ADIGATOR
PFLAG = ADIGATOR.PRINT.FLAG;
if isa(x,'cada')
  if ~isempty(x.func.value)
    str = x.func.name;
    num = x.func.value;
    ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
  else
    error('overloaded sub2ind only meant to be used on known numeric objects')
  end
else
  num = x;
  if PFLAG
    str = cadaindprint(x);
  else
    str = [];
  end
end

end