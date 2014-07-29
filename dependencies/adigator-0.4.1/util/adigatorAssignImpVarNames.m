function adigatorAssignImpVarNames(VarID,VarStr,SubsFlag)
%function adigatorAssignImpVarNames(VarID,VarStr,SubsFlag)
% this module assigns values to the ADIGATOR.VARINFO.NAMES and
% ADIGATOR.VARINFO.NAMELOCS fields for a given user defined Assignment
% Variable
% --------------------- Input Information ------------------------------- %
% VarID - .id field of the overloaded variable
% VarStr - user string which variable is named in output file
% SubsFlag - binary flag, true means variable was a subs-assignment
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
cadaStrCmp = strcmp(VarStr,ADIGATOR.VARINFO.NAMES);

if ~any(cadaStrCmp)
  %--Create New VARNAMES Cell--
  VarNameLoc = length(ADIGATOR.VARINFO.NAMES)+1;
  %--Set VARNAMELOCS--
  ADIGATOR.VARINFO.NAMES{VarNameLoc,1} = VarStr;
  ADIGATOR.VARINFO.NAMELOCS(VarID,1) = VarNameLoc;
else
  %--Already has VARNAMES Cell--
  NameIndices = 1:length(ADIGATOR.VARINFO.NAMES);
  VarNameLoc = NameIndices(cadaStrCmp);
  %--Set VARNAMELOCS--
  if length(VarNameLoc) > 1
    for Vcount = VarNameLoc(2:end)
      ADIGATOR.VARINFO.NAMES{Vcount} = '-';
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(:,1)==Vcount,1) = VarNameLoc(1);
    end
    VarNameLoc = VarNameLoc(1);
  end
  ADIGATOR.VARINFO.NAMELOCS(VarID,1) = VarNameLoc;
end
if SubsFlag == 0
    ADIGATOR.VARINFO.NAMELOCS(VarID,3) = 1;
end
return
end