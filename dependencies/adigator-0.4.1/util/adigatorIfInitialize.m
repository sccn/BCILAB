function adigatorIfInitialize(IfCount,varargin)
% This transformation routine is placed prior to an IF/ELSEIF/ELSE set in
% the intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMbros  = nargin-1;
CondVars = varargin;
if ADIGATOR.OPTIONS.PREALLOCATE
  % Pre-Allocating cells/structures - have to treat everything as if it is
  % a true statement.
elseif ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  
  % Do Naming scheme
  for Bcount = 1:NUMbros
    BroVar = CondVars{Bcount};
    if isa(BroVar,'cada')
      Vcount = BroVar.id;
      ADIGATOR.VARINFO.NAMELOCS(Vcount,3) = Inf; % Dont print derivs
      Vcount = Vcount-1;
      while ~ADIGATOR.VARINFO.NAMELOCS(Vcount,1)
        % Dont print derivs of any intermediate vars leading up to this one
        ADIGATOR.VARINFO.NAMELOCS(Vcount,3) = Inf;
        Vcount = Vcount-1;
      end
    end
  end
  if ~isa(BroVar,'cada')
    ADIGATOR.IFDATA(IfCount).ELSEFLAG = 1;
  else
    ADIGATOR.IFDATA(IfCount).ELSEFLAG = 0;
  end

  ADIGATOR.IFDATA(IfCount).BROS = ...
    struct('START',...
    cell(NUMbros,1),'END',cell(NUMbros,1),'RUNFLAG',cell(NUMbros,1),...
    'PRIORDEP',cell(NUMbros,1),'INNERDEP',cell(NUMbros,1),'REMAPS',cell(NUMbros,1));
  ADIGATOR.IFDATA(IfCount).FORLOC =...
    [ADIGATOR.FORINFO.OUTERLOC ADIGATOR.FORINFO.INNERLOC];
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %

  % ---------------- Parse If/Elseif Conditionals ----------------------- %
  if ADIGATOR.EMPTYFLAG
    % We are in an empty run - set all to empty
    ADIGATOR.IFDATA(IfCount).EMPTYFLAG = 1;
    for Bcount = 1:NUMbros
      ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
    end
    TrueFlag = 0;
  else
    ADIGATOR.IFDATA(IfCount).EMPTYFLAG = 0;
    TrueFlag  = 0;
    FalseFlag = 0;
    if ADIGATOR.IFDATA(IfCount).ELSEFLAG
      % there is an ELSE statement
      for Bcount = 1:NUMbros-1
        BroVar = CondVars{Bcount};
        if TrueFlag
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
        elseif ~isempty(BroVar.func.value)
          if logical(BroVar.func.value)
            % This is TRUE
            TrueFlag = Bcount;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          else
            % This is FALSE
            FalseFlag = 1;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
          end
        else
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
        end
      end
      if TrueFlag
        ADIGATOR.IFDATA(IfCount).BROS(NUMbros).RUNFLAG = 0;
      else
        ADIGATOR.IFDATA(IfCount).BROS(NUMbros).RUNFLAG = 1;
      end
    else
      % there is not an ELSE statement
      for Bcount = 1:NUMbros
        BroVar = CondVars{Bcount};
        if TrueFlag
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
        elseif ~isempty(BroVar.func.value)
          if logical(BroVar.func.value)
            % This is TRUE
            TrueFlag = Bcount;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          else
            % This is FALSE
            FalseFlag = 1;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
          end
        else
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
        end
      end
    end
    TrueFlag = TrueFlag | FalseFlag;
    SaveIncomingVars(IfCount);
  end
  % ------------------- Take Car of Overmapping ------------------------- %
  OvermapIncomingVars(IfCount,TrueFlag);
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  
  % Parse If/Elseif Statements to make sure there isnt some goofy case
  % where user wrote an IF/ELSEIF statement that never fires.
  if ADIGATOR.EMPTYFLAG
    % We are in an empty run - set all to empty
    ADIGATOR.IFDATA(IfCount).PRINTFLAG = 0;
    ADIGATOR.IFDATA(IfCount).EMPTYFLAG = 1;
    for Bcount = 1:NUMbros
      ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
    end
  else
    ADIGATOR.IFDATA(IfCount).EMPTYFLAG = 0;
    TrueFlag = 0;
    ADIGATOR.IFDATA(IfCount).PRINTFLAG = 0;
    if ADIGATOR.IFDATA(IfCount).ELSEFLAG
      % there is an ELSE statement
      for Bcount = 1:NUMbros-1
        BroVar = CondVars{Bcount};
        if TrueFlag
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
        elseif ~isempty(BroVar.func.value)
          if logical(BroVar.func.value)
            % This is TRUE
            TrueFlag = Bcount;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          else
            % This is FALSE
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
          end
        else
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          ADIGATOR.IFDATA(IfCount).PRINTFLAG = 1;
        end
      end
      if TrueFlag
        ADIGATOR.IFDATA(IfCount).BROS(NUMbros).RUNFLAG = 0;
      else
        ADIGATOR.IFDATA(IfCount).BROS(NUMbros).RUNFLAG = 1;
      end
    else
      % there is not an ELSE statement
      for Bcount = 1:NUMbros
        BroVar = CondVars{Bcount};
        if TrueFlag
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
        elseif ~isempty(BroVar.func.value)
          if logical(BroVar.func.value)
            % This is TRUE
            TrueFlag = Bcount;
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          else
            % This is FALSE
            ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 0;
          end
        else
          ADIGATOR.IFDATA(IfCount).BROS(Bcount).RUNFLAG = 1;
          ADIGATOR.IFDATA(IfCount).PRINTFLAG = 1;
        end
      end
    end
    SaveIncomingVars(IfCount);
  end
end
end

function SaveIncomingVars(IfCount)
global ADIGATOR ADIGATORVARIABLESTORAGE
PriorDepLocs = ADIGATOR.IFDATA(IfCount).PRIORDEP.Locs;
SaveLocs = ADIGATOR.VARINFO.SAVE.IF(:,1);


nOut = length(PriorDepLocs);

if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.RUNFLAG == 2 ...
    && ADIGATOR.IFDATA(IfCount).OUTERFLAG
  % Unrolling loops, on a conditional block which gets run twice in a row,
  % on the second evaluation.
  % - Need to reset the saved variables to what they were when we first
  % grabbed them.
  Vars = ADIGATOR.IFDATA(IfCount).PRIORDEP.Vars;
  for iOut = 1:nOut
    VarLoc = PriorDepLocs(iOut);
    SaveLoc = SaveLocs(VarLoc);
    ADIGATORVARIABLESTORAGE.SAVE{SaveLoc} = Vars{iOut};
  end
else
  % When a conditional block is dependent upon a var occuring prior to
  % loop, but that var gets redefined in a branch and then used in a later
  % branch, we save all the prior occuring vars at this point so as the
  % save locs may get assigned new objects as we evaluate the block. The
  % IfIterStart routine may then access them and return them to the
  % evaluating workspace at the appropriate time.
  SaveVars = ADIGATORVARIABLESTORAGE.SAVE;
  Vars = cell(1,nOut);
  for iOut = 1:nOut
    VarLoc = PriorDepLocs(iOut);
    SaveLoc = SaveLocs(VarLoc);
    Vars{iOut} = SaveVars{SaveLoc};
  end
  ADIGATOR.IFDATA(IfCount).PRIORDEP.Vars = Vars;
end

end

function OvermapIncomingVars(IfCount,TrueFlag)
global ADIGATOR ADIGATORVARIABLESTORAGE
if TrueFlag
  % One of the statements returned true or false, dont necessarily have to
  % do all of the Prior Overmaps
  nBros = length(ADIGATOR.IFDATA(IfCount).BROS);
  RunFlags = zeros(nBros,1);
  for BroCount = 1:nBros
    if ADIGATOR.IFDATA(IfCount).BROS(BroCount).RUNFLAG
      RunFlags(BroCount) = BroCount;
    end
  end
  BroCounts = nonzeros(RunFlags);
  nRuns     = length(BroCounts);
  if ~nRuns; TrueFlag = 0; end;
  RunCounts = cell(nRuns,1);
  for rCount = 1:nRuns
    BroCount = BroCounts(rCount);
    Start = ADIGATOR.IFDATA(IfCount).BROS(BroCount).START;
    End   = ADIGATOR.IFDATA(IfCount).BROS(BroCount).END;
    RunCounts{rCount} = Start:End;
  end
end
Overmap = ADIGATOR.IFDATA(IfCount).OVERMAP;
for Ocount = 1:length(Overmap)
  PriorCount = Overmap(Ocount).Prior;
  if ~isempty(PriorCount)
    if TrueFlag
      % See if we need to overmap
      NameLoc = ADIGATOR.VARINFO.NAMELOCS(PriorCount,1);
      OverFlag = 0;
      for rCount = 1:nRuns
        if ~any(ADIGATOR.VARINFO.NAMELOCS(RunCounts{rCount},1) == NameLoc)
          OverFlag = 1;
          break
        end
      end
      if ~OverFlag
        continue
      end
    end
    % ------------------- Overmap Incoming Variables -------------------- %
    OverLoc = Overmap(Ocount).Outer(1);
    SaveLoc = ADIGATOR.VARINFO.SAVE.IF(PriorCount,1);
    SaveVar = ADIGATORVARIABLESTORAGE.SAVE{SaveLoc};
    ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc} = SaveVar;
  elseif ADIGATOR.FORINFO.FLAG
    % Remove the overmap (in case of for loop)
    OverLoc = Overmap(Ocount).Outer(1);
    ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc} = [];
  end
end
end