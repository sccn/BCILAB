function [outEvalStr, outEvalVar] = adigatorIfIterStart(IfCount,BroCount)
% This transformation routine is placed prior the the beginning of each
% IF/ELSEIF/ELSE block in the intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.OPTIONS.PREALLOCATE
  % Pre-Allocating cells/structures - have to treat everything as if it is
  % a true statement.
  outEvalStr = [];
  outEvalVar = [];
elseif ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  ADIGATOR.IFDATA(IfCount).BROS(BroCount).START = ADIGATOR.VARINFO.COUNT;
  outEvalVar = [];
  outEvalStr = [];
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  ADIGATOR.EMPTYFLAG = ~ADIGATOR.IFDATA(IfCount).BROS(BroCount).RUNFLAG;
  if ~ADIGATOR.EMPTYFLAG  && ~isempty(ADIGATOR.IFDATA(IfCount).BROS(BroCount).PRIORDEP)
    [outEvalStr outEvalVar] = GetIncomingVariables(...
      ADIGATOR.IFDATA(IfCount).PRIORDEP,ADIGATOR.IFDATA(IfCount).BROS(BroCount).PRIORDEP);
  else
    outEvalVar = [];
    outEvalStr = [];
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  if BroCount == 1
    % ---------------------------- IF ----------------------------------- %
    if ADIGATOR.IFDATA(IfCount).BROS(1).RUNFLAG && ADIGATOR.IFDATA(IfCount).PRINTFLAG
      fprintf(fid,[indent,'if cadaconditional1\n']);
    end
  else
    % ------------------------ ELSEIF/ELSE ------------------------------ %
    if ADIGATOR.IFDATA(IfCount).BROS(BroCount).RUNFLAG
      if ADIGATOR.EMPTYFLAG
        % Previous Case didnt run, check to make sure IF statement has been
        % printed.
        IfFlag = 0;
        for Bcount = 1:BroCount-1
          if ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG
            IfFlag = 1; break
          end
        end
      else
        IfFlag = 1;
      end
      if IfFlag
        if ~ADIGATOR.IFDATA(IfCount).ELSEFLAG ||...
            BroCount < length(ADIGATOR.IFDATA(IfCount).BROS)
          fprintf(fid,[indent,'elseif cadaconditional%1.0d\n'],BroCount);
        else
          fprintf(fid,[indent,'else\n']);
        end
      else
        if (~ADIGATOR.IFDATA(IfCount).ELSEFLAG ||...
            BroCount < length(ADIGATOR.IFDATA(IfCount).BROS)) ...
            && ADIGATOR.IFDATA(IfCount).PRINTFLAG
          fprintf(fid,[indent,'if cadaconditional%1.0d\n'],BroCount);
        end
      end
    end
  end
  ADIGATOR.PRINT.INDENT = [indent,'    '];
  ADIGATOR.EMPTYFLAG = ~ADIGATOR.IFDATA(IfCount).BROS(BroCount).RUNFLAG;
  if ~ADIGATOR.EMPTYFLAG  &&  ~isempty(ADIGATOR.IFDATA(IfCount).BROS(BroCount).PRIORDEP)
    [outEvalStr outEvalVar] = GetIncomingVariables(...
      ADIGATOR.IFDATA(IfCount).PRIORDEP,ADIGATOR.IFDATA(IfCount).BROS(BroCount).PRIORDEP);
  else
    outEvalVar = [];
    outEvalStr = [];
  end
end

end

function [outEvalStr outEvalVar] = GetIncomingVariables(PriorDepStruc,PriorDepLocs)
global ADIGATOR

NameLocs = ADIGATOR.VARINFO.NAMELOCS(:,1);
Names    = ADIGATOR.VARINFO.NAMES;

nOut = length(PriorDepLocs);
outEvalStr = cell(nOut,1);
outEvalVar = cell(nOut,1);
for iOut = 1:nOut
  VarLoc = PriorDepLocs(iOut);
  NameLoc = NameLocs(VarLoc);
  VarName = Names{NameLoc};
  outEvalStr{iOut} = sprintf([VarName,' = adigatorIfEvalVar{%1.0f};'],iOut);
  outEvalVar{iOut} = PriorDepStruc.Vars{PriorDepStruc.Locs == VarLoc};
end
end
