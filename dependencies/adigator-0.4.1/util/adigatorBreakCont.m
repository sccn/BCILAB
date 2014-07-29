function adigatorBreakCont(type,IfCount,BroCount)
% All user defined BREAKS/CONTINUES are replaced with this function in the
% intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

if ~ADIGATOR.RUNFLAG
  if ~ADIGATOR.FORINFO.FLAG
    error('breaks and continues must be contained within a loop')
  elseif ADIGATOR.OPTIONS.UNROLL
    error('breaks and continues cannot be used when unrolling loops.')
  end
  if strcmp(type,'break')
    ADIGATOR.BREAKLOCS(end+1,:) = [IfCount, BroCount,...
      ADIGATOR.FORINFO.INNERLOC,ADIGATOR.FORINFO.OUTERLOC];
  else
    ADIGATOR.CONTLOCS(end+1,:)  = [IfCount, BroCount,...
      ADIGATOR.FORINFO.INNERLOC,ADIGATOR.FORINFO.OUTERLOC];
  end
elseif ADIGATOR.RUNFLAG == 2
  fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,type,'\n']);
end