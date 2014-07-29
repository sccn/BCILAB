function LoopVari = adigatorForIterStart(ForCount,ForIter)
% This transformation routine is called at the start of each FOR loop
% iteration in the intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFORDATA
if ADIGATOR.OPTIONS.PREALLOCATE
  LoopVari = ForIter;
  return
end
NUMvod = ADIGATOR.NVAROFDIFF;
if ~ADIGATOR.OPTIONS.UNROLL
  ADIGATORFORDATA(ForCount).COUNT.SUBSREF    = 0;
  ADIGATORFORDATA(ForCount).COUNT.SUBSASGN   = 0;
  ADIGATORFORDATA(ForCount).COUNT.SPARSE     = 0;
  ADIGATORFORDATA(ForCount).COUNT.NONZEROS   = 0;
  ADIGATORFORDATA(ForCount).COUNT.HORZCAT    = 0;
  ADIGATORFORDATA(ForCount).COUNT.VERTCAT    = 0;
  ADIGATORFORDATA(ForCount).COUNT.TRANSPOSE  = 0;
  ADIGATORFORDATA(ForCount).COUNT.RESHAPE    = 0;
  ADIGATORFORDATA(ForCount).COUNT.REPMAT     = 0;
  ADIGATORFORDATA(ForCount).COUNT.SIZE       = 0;
  ADIGATORFORDATA(ForCount).COUNT.OTHER      = 0;
end
ADIGATORFORDATA(ForCount).COUNT.ITERATION  = ForIter;

func      = struct('name',[],'size',[1 1],'zerolocs',[],'value',ForIter);
func.name = ADIGATORFORDATA(ForCount).COUNTNAME;
deriv     = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.RUNFLAG
  ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).START;
  ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  LoopVari = ForIter;
  return
elseif ADIGATOR.RUNFLAG == 1
  ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).START;
  ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  if ~isempty(ADIGATOR.VARINFO.OVERMAP.CONT)
    % At the start of each iteration we need to clear any CONTINUE OVERMAPS
    % associated with this loop.
    Start = ADIGATORFORDATA(ForCount).START;
    End   = ADIGATORFORDATA(ForCount).END;
    ContOverLocs = nonzeros(unique(ADIGATOR.VARINFO.OVERMAP.CONT(Start:End)));
    if ADIGATORFORDATA(ForCount).PARENTLOC && ...
        (ADIGATORFORDATA(ForCount).PARENTLOC > 1 || ~ADIGATOR.FORINFO.FUNLOOP)...
        && ForIter == 1
      % If this is a nested loop and a the parent loop has a variable prior
      % to this loop which shares a CONTINUE OVERMAP with variables in this
      % loop, then we need to 
      % 1. If on first itertion, Save the current CONTINUE OVERMAP
      % 2. Clear the CONTINUE OVERMAP
      % 3. adigatorForIterEnd will take care of making the CONTINUE overmap
      % proper again.
      PrevCounts      = 1:Start-1;
      ContRestoreLocs = ContOverLocs;
      for COcount = 1:length(ContOverLocs)
        COloc = ContOverLocs(COcount);
        if ~any(ADIGATOR.VARINFO.OVERMAP.CONT(PrevCounts) == COloc)
          ContRestoreLocs(COcount) = 0;
        end
      end
      ContOverLocs    = nonzeros(ContOverLocs);
      ContRestoreLocs = nonzeros(ContRestoreLocs);
    else
      ContRestoreLocs = [];
    end
    if isempty(ContRestoreLocs)
      ClearOvermaps(ContOverLocs);
    else
      RestoreOvermaps = ClearOvermaps(ContOverLocs,ContRestoreLocs);
      ADIGATORFORDATA(ForCount).CONTRESTORE.LOCS = ContRestoreLocs;
      ADIGATORFORDATA(ForCount).CONTRESTORE.VARS = RestoreOvermaps;
    end
  end
  if ~isempty(ADIGATOR.VARINFO.OVERMAP.BREAK) && ForIter == 1
    % Need to clear all of the BREAK overmaps for this loop
    Start = ADIGATORFORDATA(ForCount).START;
    End   = ADIGATORFORDATA(ForCount).END;
    BreakOverLocs = nonzeros(unique(abs(ADIGATOR.VARINFO.OVERMAP.BREAK(Start:End))));
    if ADIGATORFORDATA(ForCount).PARENTLOC && ...
        (ADIGATORFORDATA(ForCount).PARENTLOC > 1 || ~ADIGATOR.FORINFO.FUNLOOP)
      % If this is a nested loop and the parent loop has a variable prior
      % to this loop which shares a BREAK OVERMAP with variables in this
      % loop, then we just dont clear the BREAK OVERMAP
      PrevCounts = 1:Start-1;
      for BOcount = 1:length(BreakOverLocs)
        BOloc = BreakOverLocs(BOcount);
        if any(abs(ADIGATOR.VARINFO.OVERMAP.BREAK(PrevCounts)) == BOloc)
          BreakOverLocs(BOcount) = 0;
        end
      end
      BreakOverLocs = nonzeros(BreakOverLocs);
    end
    ClearOvermaps(BreakOverLocs);
  end
elseif ADIGATOR.RUNFLAG == 2
  ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).START;
  func.value = [];
end
LoopVari  = cada([],func,deriv);

end

function varargout = ClearOvermaps(OverLocs,varargin)
global ADIGATORVARIABLESTORAGE

if nargin == 2
  RestoreLocs = varargin{1};
  varargout{1} = ADIGATORVARIABLESTORAGE.OVERMAP(RestoreLocs);
end
ADIGATORVARIABLESTORAGE.OVERMAP(OverLocs) = cell(1);

end