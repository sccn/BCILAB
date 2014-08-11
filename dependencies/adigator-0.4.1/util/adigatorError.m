function adigatorError(IfCount,BroCount,ErrorMsg)
% All user error functions are replaced with this function in the
% intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

if ~ADIGATOR.RUNFLAG
  ADIGATOR.ERRORLOCS(end+1,:) = [IfCount, BroCount];
elseif ADIGATOR.RUNFLAG == 2 && ~ADIGATOR.EMPTYFLAG
  fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,ErrorMsg,'\n']);
end