function ifLoop = adigatorIfLooper(IfCount)
% This transformation routine is called when we are unrolling loops and
% encounter an IF statement.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.OPTIONS.PREALLOCATE
  ifLoop = 1;
elseif ~ADIGATOR.RUNFLAG
  ifLoop = 0;
  ADIGATOR.IFDATA(IfCount).OUTERFLAG = 1;
elseif ~ADIGATOR.IFINFO.INNERLOC
  % Outermost IF statement
  ADIGATOR.IFDATA(IfCount).RESETLOC = 0;
  ADIGATOR.IFINFO.INNERLOC = IfCount;
  ifLoop = [1 2];
elseif ADIGATOR.IFINFO.INNERLOC
  ADIGATOR.IFDATA(IfCount).RESETLOC = ADIGATOR.IFINFO.INNERLOC;
  ADIGATOR.IFINFO.INNERLOC = IfCount;
  if ADIGATOR.RUNFLAG == 1
    % Overmapping evaluation of outer conditional
    ifLoop = 1;
  else
    % Printing an outer conditional statement, but there exists a loop in
    % between the outer conditional and this conditional.
    ifLoop = [1 2];
  end
end