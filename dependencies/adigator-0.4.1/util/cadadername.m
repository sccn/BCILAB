function derivstr = cadadername(funcname,Vcount,varargin)
% This gets the derivative name of some overloaded object.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if nargin == 2
    VarID = ADIGATOR.VARINFO.COUNT;
elseif nargin == 3
    VarID = varargin{1};
end

if ADIGATOR.PRINT.FLAG == 1
  nameindex = ADIGATOR.VARINFO.NAMELOCS(VarID);
  if nameindex && (ADIGATOR.VARINFO.NAMELOCS(VarID,2)>=0)
    if ADIGATOR.DERNUMBER == 1
      derivstr = [ADIGATOR.VARINFO.NAMES{nameindex},'.d',...
        ADIGATOR.VAROFDIFF(Vcount).name];
    elseif length(funcname) > 2 && strcmp(funcname(end-1:end),'.f')
      derivstr = [funcname(1:end-2),'.d',ADIGATOR.VAROFDIFF(Vcount).name];
    else
      derivstr = [funcname,'d',ADIGATOR.VAROFDIFF(Vcount).name];
    end
  elseif length(funcname) > 2 && strcmp(funcname(end-1:end),'.f')
    derivstr = [funcname(1:end-2),'.d',ADIGATOR.VAROFDIFF(Vcount).name];
  else
    derivstr = [funcname,'d',ADIGATOR.VAROFDIFF(Vcount).name];
  end
else
  derivstr = [funcname,'d',ADIGATOR.VAROFDIFF(Vcount).name];
end
