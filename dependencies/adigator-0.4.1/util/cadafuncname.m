function [funcstr,derivflag] = cadafuncname(varargin)
% Get function string of an object and a flag stating whether or not we
% need to calculate derivatives.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

if nargin == 0
  VarID = ADIGATOR.VARINFO.COUNT;
elseif nargin == 1
  VarID = varargin{1};
end
derivflag = 0;
if ADIGATOR.PRINT.FLAG

  nameindex = ADIGATOR.VARINFO.NAMELOCS(VarID,1);

  if nameindex == 0
    % intermediate variable
    namenum = ADIGATOR.VARINFO.NAMELOCS(VarID,2);
    funcstr = sprintf('cada%1.0df%1.0d',ADIGATOR.DERNUMBER,namenum);
    if ~isinf(ADIGATOR.VARINFO.NAMELOCS(VarID,3))
      derivflag = 1;
    end
  else
    % important variable
    if ADIGATOR.DERNUMBER == 1
      if isinf(ADIGATOR.VARINFO.NAMELOCS(VarID,3)) || ...
          any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(VarID,1)...
          ==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
        funcstr = ADIGATOR.VARINFO.NAMES{nameindex};
        derivflag = 0;
      elseif ADIGATOR.VARINFO.NAMELOCS(VarID,2) < 0
        funcstr = ADIGATOR.VARINFO.NAMES{nameindex};
        CheckStrings = ADIGATOR.DERIVCHECKS(-ADIGATOR.VARINFO.NAMELOCS(VarID,2)).STRINGS;
        if isempty(CheckStrings) || any(~(cellfun('isempty',regexp(funcstr,CheckStrings,'start'))))
          derivflag = 1;
        end
      else
        funcstr = [ADIGATOR.VARINFO.NAMES{nameindex},'.f'];
        derivflag = 1;
      end
    else
      funcstr = [ADIGATOR.VARINFO.NAMES{nameindex}];
      if ADIGATOR.DERIVCHECKS(ADIGATOR.FILE.FUNID).NUM
        % There are only certain variables that we need to print out
        % derivatives for.
        CheckStrings = ADIGATOR.DERIVCHECKS(ADIGATOR.FILE.FUNID).STRINGS;
        % Check for variables ending in CheckStrings
        if isempty(CheckStrings) || any(~(cellfun('isempty',regexp(funcstr,CheckStrings,'start'))))
          derivflag = 1;
        end
      else
        derivflag = 1;
      end

    end
  end
else
  funcstr = sprintf('f%1.0d',VarID);
  derivflag = 0;
end

end