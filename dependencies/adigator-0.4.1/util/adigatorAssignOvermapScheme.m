function [FORDATA, IFDATA, VARINFO,VARSTORAGE] = adigatorAssignOvermapScheme(FunID,...
    FunAsLoopFlag,FORDATA,IFDATA,VARINFO,BREAKLOCS,CONTLOCS,ERRORLOCS,UNROLL)
% This function is used to Assign the Overmap scheme for all FOR loops and
% IF/ELSEIF/ELSE statements. All of the details on the overmapping scheme
% may be found in the comments of this code.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
NAMELOCS = VARINFO.NAMELOCS;
LASTOCC  = VARINFO.LASTOCC;
nTotalVars = VARINFO.COUNT-1;

if ~isempty(IFDATA)
  OVERMAP.IF  = zeros(nTotalVars,5);
  SAVE.IF     = zeros(nTotalVars,2);
  RETURN      = zeros(nTotalVars,1);
  LOOPMAP     = zeros(nTotalVars,1);
  if ~isempty(FORDATA)
    for Lcount = 1:length(FORDATA)
      Lstart = FORDATA(Lcount).START;
      LOOPMAP(Lstart) = Lcount;
    end
  end
else
  OVERMAP.IF  = [];
  SAVE.IF     = [];
  RETURN      = [];
end
if ~isempty(FORDATA) && ~UNROLL && (FunID == 1 || FunAsLoopFlag || length(FORDATA) > 1)
  OVERMAP.FOR = zeros(nTotalVars,2);
  SAVE.FOR    = zeros(nTotalVars,2);
else
  OVERMAP.FOR = [];
  SAVE.FOR    = [];
end
if ~isempty(CONTLOCS)
  OVERMAP.CONT = zeros(nTotalVars,1);
else
  OVERMAP.CONT = [];
end
if ~isempty(BREAKLOCS) || ~isempty(CONTLOCS)
  OVERMAP.BREAK = zeros(nTotalVars,1);
else
  OVERMAP.BREAK = [];
end
SaveCount = 0;
OverCount = 0;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ LOOP STATEMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~ %%                 
if ~UNROLL
  % if OVERMAP.FOR(i,1) = j, then variable i is contained within a FOR loop and
  %   shares the overmap variable contained in
  %   ADIGATORVARIABLESTORAGE.OVERMAP{j}
  % if OVERMAP.FOR(i,2) = j, then variable i is outside of a FOR loop and the
  %   gets used within the FOR loop, then changed after an iteration, the
  %   overmap for this is then contained in slot j.
  % if SAVE.FOR(i,1) = j, then variable i is active previous to a FOR loop
  %   is redefined within the loop. In this case we store the variable, then
  %   once entering the FOR loop, we union it with the proper variable from
  %   within the loop.
  % if SAVE.FOR(i,2) = j, then variable i is defined within a FOR loop and is
  %   used outside of the FOR loop - we save this so we know what to Re-Map
  %   to on the printing run.

  for ForCount = 1:length(FORDATA)
    if ForCount == 1 && FunID > 1 && ~FunAsLoopFlag
      continue
    end
    Start = FORDATA(ForCount).START;
    End   = FORDATA(ForCount).END;
    AllVarCounts  = Start:End;
    AllAsgnCounts = AllVarCounts(logical(NAMELOCS(AllVarCounts,1)));
    
    % -------------------- Outer Loop - Do Assignments -------------------- %
    if ~FORDATA(ForCount).PARENTLOC
      % Is an outer loop - Define all of the overmapping within.
      OuterStart = Start;
      for Vcount = 1:length(AllAsgnCounts)
        AAcount  = AllAsgnCounts(Vcount);
        SubsFlag = NAMELOCS(AAcount,3);
        NameLoc  = NAMELOCS(AAcount,1);
        if ~SubsFlag
          % Is a SubsAsgn - Check to See if the Corresponding Direct
          % Assignment has an OverMap
          PrevAsgnCounts = AllAsgnCounts(1:Vcount-1);
          DAcount = PrevAsgnCounts(find(NAMELOCS(PrevAsgnCounts,1) == NameLoc & ...
            NAMELOCS(PrevAsgnCounts,3),1,'last'));
          % DAcount is the last direct assignment within the loop that
          % shares the same name.
          if ~isempty(DAcount)
            % There is a direct assignment within the loop
            OVERMAP.FOR(AAcount,1) = OVERMAP.FOR(DAcount,1);
          else
            % There is not a direct assignment within the loop which shares
            % the same name as this variable
            SAcount = PrevAsgnCounts(find(NAMELOCS(PrevAsgnCounts,1) == NameLoc,1,'last'));
            if ~isempty(SAcount)
              % There was a subs-assignment to same variable prior to this
              % within the loop, use that overmap.
              OVERMAP.FOR(AAcount,1) = OVERMAP.FOR(SAcount,1);
            else
              % Give this a new overmap location
              OverCount = OverCount+1;
              OVERMAP.FOR(AAcount,1) = OverCount;
            end
          end
        else
          % Direct Assignment - Give this its own OverMap loc
          OverCount = OverCount+1;
          OVERMAP.FOR(AAcount) = OverCount;
        end
      end
      
      % -------------------- Set Variables to Save ------------------------ %
      % Save all objects created within the loop which get used after the
      % loop
      if ~FunAsLoopFlag
        SaveCounts = AllAsgnCounts(LASTOCC(AllAsgnCounts,1) > AllAsgnCounts(end));
        for Scount = 1:length(SaveCounts)
          SaveLoc = SAVE.FOR(SaveCounts(Scount),1);
          if ~SaveLoc
            SaveCount = SaveCount+1;
            SaveLoc   = SaveCount;
          end
          SAVE.FOR(SaveCounts(Scount),2) = SaveLoc;
        end
      end
      % ------- Overmap Starting Variables with Ending Variables ---------- %
      % (This is outer loop, in the main function the starting variables get
      % overmapped when they come into the outer loop)
      if ~FunAsLoopFlag
        NameLocs      = unique(NAMELOCS(AllAsgnCounts,1));
        PrevVarCounts = 1:Start-1;
        for Ncount = 1:length(NameLocs)
          NameLoc   = NameLocs(Ncount);
          VarCount  = PrevVarCounts(find(NAMELOCS(PrevVarCounts,1) == NameLoc,1,'last'));
          if ~isempty(VarCount) && LASTOCC(VarCount,1) >= Start
            SaveCount = SaveCount+1;
            SAVE.FOR(VarCount,1) = SaveCount;
            FORDATA(ForCount).PREVOVERMAP =...
              [FORDATA(ForCount).PREVOVERMAP VarCount];
            OverLoc = OVERMAP.FOR(AllAsgnCounts(find(NAMELOCS(AllAsgnCounts,1)...
              ==NameLoc,1,'last')),1);
            OVERMAP.FOR(VarCount,2) = OverLoc;
          end
        end
      end
    else
      % ------- Overmap Starting Variables with Ending Variables ---------- %
      % (This is some nested loop, make it such that all variables coming in
      % share an overmap with variables at the end of the loop - will make
      % some empty overmap spots but thats okay)
      NameLocs      = unique(NAMELOCS(AllAsgnCounts,1));
      PrevVarCounts = 1:Start-1;
      for Ncount = 1:length(NameLocs)
        NameLoc   = NameLocs(Ncount);
        VarCount  = PrevVarCounts(find(NAMELOCS(PrevVarCounts,1) == NameLoc,1,'last'));
        if ~isempty(VarCount) && LASTOCC(VarCount,1) >= Start && VarCount >= OuterStart
          SaveCount = SaveCount+1;
          SAVE.FOR(VarCount,1) = SaveCount;
          OverLoc = OVERMAP.FOR(AllAsgnCounts(find(NAMELOCS(AllAsgnCounts,1)...
            ==NameLoc,1,'last')),1);
          OldOverLoc = OVERMAP.FOR(VarCount,1);
          OVERMAP.FOR(OVERMAP.FOR == OldOverLoc) = OverLoc;
        end
      end
    end
  end
end
  
%% ~~~~~~~~~~~~~~~~~~~~~ CONDITIONAL STATEMENTS ~~~~~~~~~~~~~~~~~~~~~~~~ %%
% ---------------------------- VARINFO STUFF ---------------------------- %
% -if OVERMAP.IF(i,1) = j, then variable i must be unioned to create overmap
% in j, OVERMAP.IF(i,2) and OVERMAP.IF(i,3) will also be built such that
% OVERMAP.IF(i,2) = k, OVERMAP.IF(i,3) = l, where the overmap j is for 
% IF block number k, branch number l
% -if OVERMAP.IF(i,4) = j, then the object i belongs to overmap j, but is
% defined prior to the set.
% -if OVERMAP.IF(i,5) = j, then the loop number j is nested within the
% branch OVERMAP.IF(i,3), and OVERMAP.IF(i,1) is not zero.
%
% This gets changed to OVERMAP.IF = OVERMAP.IF(:,[1 5]) upon exit of this
% subfunction, the columns 2:4 are only used for building the overmapping
% scheme.
%
% -if SAVE.IF(i,1) = j, then variable occurs prior to the conditional set
% and is either saved for the overmap or because the variable is changed
% within one branch and another branch is dependent upon it.
% -if SAVE.IF(i,2) = j, then this variable is needed for an overmap, but it
% gets used within its branch after the last assignment to it. Save this in
% the PRINTING evaluation so that the brother IfIterEnd may do the remap.
%
% -if RETURN(i) = 1, then return variable i as its overmap.
%
% --------------------------- IFDATA STUFF ------------------------------ %
% 
% IFDATA(i).OVERMAP(j).Prior - the count of the variable occuring prior
% to block but belonging to overmap j.
%
% IFDATA(i).OVERMAP(j).Outer = [o s]
%     o: location of overmap
%     k: if k ~= 0, then must save overmap in save location k (for a
%     parent) in the printing evaluation
%
% IFDATA(i).OVERMAP(j).Inner = array of, oi, where oi is the location of a
% child overmap which belongs to this one. Need to bring this into the
% overmap during the overmapping eval.
%
% IFDATA(i).OVERMAP(j).Return - if not empty, then need to return this
% overmap to the evaluating workspace from IfIterEnd (in both evaluations).
% The value stored here will be the NameLoc of the variable which we assign
% the overmap to.
%
% IFDATA(i).BROS(b).PRIORDEP = array of variable counts, cj, where each
% variable must be loaded prior to the evaluation of branch b. Each
% variable will be saved in the location SAVE.IF(cj,1). Used in both
% overmapping and printing evaluation.
% 
% IFDATA(i).BROS(b).REMAPS = 2 x n array of integers, for column j=[j1;j2]
%       if j2 = 0, then must remap object saved in location
%       IFDATA(i).OVERMAP(cj).Prior to the overmap in location
%       IFDATA(i).OVERMAP(cj).Outer(1).
%       else must remap object saved in location j1 to the overmap in
%       location OVERMAP.IF(j2)

for IFcount = 1:length(IFDATA)
  BroData = IFDATA(IFcount).BROS;
  nBros   = length(BroData);
  Start   = BroData(1).START;
  End     = BroData(nBros).END;
  AllVarCounts  = Start:End;
  ALLVarCounts  = Start:End;
  OUTERLOC = IFDATA(IFcount).FORLOC(1);
  INNERLOC = IFDATA(IFcount).FORLOC(2);
  PRIORDEP = IFDATA(IFcount).BROS(1).PRIORDEP;
  for BroCount = 2:nBros
    PRIORDEP = [PRIORDEP IFDATA(IFcount).BROS(BroCount).PRIORDEP]; %#ok<AGROW>
  end
  % --------- Remove Any Branches with BREAK/CONTINUE/ERROR ------------- %
  BreakLocs = BREAKLOCS(BREAKLOCS(:,1) == IFcount,:);
  for BreakCount = 1:size(BreakLocs,1)
    BroBreakLoc = BreakLocs(BreakCount,2);
    bStart      = BroData(BroBreakLoc).START;
    bEnd        = BroData(BroBreakLoc).END;
    AllVarCounts(bStart:bEnd) = 0;
  end
  ContLocs  = CONTLOCS(CONTLOCS(:,1) == IFcount,:);
  for ContCount = 1:size(ContLocs,1)
    BroContLoc = ContLocs(ContCount,2);
    bStart     = BroData(BroContLoc).START;
    bEnd       = BroData(BroContLoc).END;
    AllVarCounts(bStart:bEnd) = 0;
  end
  ErrorLocs  = ERRORLOCS(ERRORLOCS(:,1) == IFcount,:);
  for ErrCount = 1:size(ErrorLocs,1)
    BroErrLoc  = ErrorLocs(ErrCount,2);
    bStart     = BroData(BroErrLoc).START;
    bEnd       = BroData(BroErrLoc).END;
    AllVarCounts(bStart:bEnd) = 0;
  end
  if ~isempty(BreakLocs) || ~isempty(ContLocs) || ~isempty(BreakLocs)
    AllVarCounts = nonzeros(AllVarCounts).';
  end
  
  % ---------------- Find all Vars we need to Overmap ------------------- %
  if ~OUTERLOC || (FunID > 1 && INNERLOC == 1)
    % All Variables which get used after the conditional set
    UsedVarCounts = AllVarCounts(LASTOCC(AllVarCounts,1) > End);
    % Name Locs of all variables assigned within the conditional set that get
    % used after the set
    AsgnVarNameLocs = unique(NAMELOCS(UsedVarCounts,1));
  else
    % We cant just check what gets used after the IF statement, since we
    % are in a FOR loop, the IF statement may redefine something which is
    % used within the FOR loop prior to the conditional set.
    AllVarNameLocs  = unique(nonzeros(NAMELOCS(AllVarCounts,1)));
    % Loop on these to see if we need to do an overmap for them. Variables
    % which get used after the conditional set need to be overmapped as
    % well as any variables which are used within the FOR loop. We can
    % check to see if a variable is used within the FOR loop by finding the
    % last occurence of the variable of the same name prior to the FOR
    % loop.
    
    % See what loops this conditional set is embedded within
    iForLoc  = INNERLOC;
    ForLoops = zeros(INNERLOC+1-OUTERLOC,1);
    ForLoops(1) = iForLoc;
    for ifCount = 2:INNERLOC+1-OUTERLOC
      iForLoc = FORDATA(iForLoc).PARENTLOC;
      if ~isempty(iForLoc)
        ForLoops(ifCount) = iForLoc;
      else
        ForLoops = ForLoops(1:ifCount-1);
        break
      end
    end
    if FunID > 1 && FunAsLoopFlag; ForLoops(end) = []; end

    for Vcount = 1:length(AllVarNameLocs)
      NameLoc = AllVarNameLocs(Vcount);
      if any(LASTOCC(AllVarCounts(NAMELOCS(AllVarCounts,1) == NameLoc)) > End)
        % Gets used after the statement - use it.
        continue
      else
        % We need to check within the FOR loops this is nested in, for a
        % given loop
        % IF
        % 1. The variable exists prior to the loop
        % 2. The variable that exists prior to the loop gets used within
        %    the loop
        % 3. The variable does not get reassigned within the loop after the
        %    conditional set.
        % THEN we need to use the variable
        NeedFlag = 0;
        for ForLoc = ForLoops
          ForStart   = FORDATA(ForLoc).START;
          ForEnd     = FORDATA(ForLoc).END;
          PrevCounts = 1:ForStart-1;
          PrevCount  = PrevCounts(find(NAMELOCS(PrevCounts,1)==NameLoc,1,'last'));
          if ~isempty(PrevCount) && LASTOCC(PrevCount) > ForStart
            % meet first two conditions, check third
            AfterCounts = End+1:ForEnd;
            if ~any(NAMELOCS(AfterCounts,1)==NameLoc)
              NeedFlag = 1;
              break
            end
          end
        end
        if ~NeedFlag
          AllVarNameLocs(Vcount) = 0;
        end
      end
    end
    AsgnVarNameLocs = nonzeros(AllVarNameLocs);
  end
  % -------------------- Set Overmap and Save Locations  ---------------- %
  nAsgnVarnameLocs = length(AsgnVarNameLocs);
  if ~INNERLOC || UNROLL
    RemapFlag  = 1;
    BROREMAPS  = cell(nBros,1);
    for BroCount = 1:nBros
      BROREMAPS{BroCount} = zeros(2,nAsgnVarnameLocs);
    end
  else
    RemapFlag  = 0;
  end
  IfOvermap = struct('Prior',cell(nAsgnVarnameLocs,1),...
    'Inner',cell(nAsgnVarnameLocs,1),'Outer',cell(nAsgnVarnameLocs,1),...
    'Return',cell(nAsgnVarnameLocs,1));
  for Vcount = 1:nAsgnVarnameLocs
    OverCount   = OverCount+1;
    OverLoc     = OverCount;
    NameLoc     = AsgnVarNameLocs(Vcount);
    ReturnCount = ALLVarCounts(find(NAMELOCS(ALLVarCounts,1) == NameLoc,1,'last'));
    IfOvermap(Vcount).Outer = [OverLoc 0];
    if ~any(AllVarCounts==ReturnCount)
      % The last assignment was made within a break/continue/error branch,
      % return this at last IfIterEnd
      IfOvermap(Vcount).Return   = NameLoc;
    end
    PrevFlag  = 0;
    for BroCount = 1:nBros
      bStart    = BroData(BroCount).START;
      bEnd      = BroData(BroCount).END;
      BroCounts = bStart:bEnd;
      LastOp    = BroCounts(find(NAMELOCS(BroCounts,1) == NameLoc,1,'last'));
      OuterOverLoc = 0;
      % LastOp = Last variable assigned within this branch
      if ~isempty(LastOp)
        % --------------------  Set IF Overmap Location ----------------- %
        if OVERMAP.IF(LastOp,1)
          % LastOp belongs to a parent conditional overmap
          OuterOverLoc            = OVERMAP.IF(LastOp,1);
          OuterIfCount            = OVERMAP.IF(LastOp,2);
          OuterOverCount          = OVERMAP.IF(LastOp,3);
          OuterSaveLoc            = SAVE.IF(LastOp,2);
          % Set this (child one)
          IfOvermap(Vcount).Outer = [OverLoc OuterSaveLoc];
          % Give info to parent
          IFDATA(OuterIfCount).OVERMAP(OuterOverCount).Inner = ...
            [IFDATA(OuterIfCount).OVERMAP(OuterOverCount).Inner OverLoc];
          % Set OVERMAP and SAVE
          OVERMAP.IF(LastOp,1:3)    = [OverLoc IFcount Vcount];
          OVERMAP.IF(LastOp,5)      = 0;
          SAVE.IF(LastOp,2)       = 0;
          LastIfCount = LastOp;
        else
          OVERMAP.IF(LastOp,1:3)    = [OverLoc IFcount Vcount];
          LastIfCount = LastOp;
        end
        
        % --------------- Check for LOOP inside of this branch ---------- %
        if any(LOOPMAP(BroCounts))
          % There is a loop within this branch - set it so that the overmap
          % is only built on the last iteration of all the loops nested
          % within.
          LoopCounts = nonzeros(LOOPMAP(BroCounts));
          OuterLoop = LoopCounts(1);
          if LastOp >= FORDATA(OuterLoop).START && LastOp <= FORDATA(OuterLoop).END
            OVERMAP.IF(LastOp,5) = OuterLoop;
          end
        end
        
        % -- Set Where Overmap Gets Returned and Remap Gets Performed --- %
        if RETURN(LastOp,1) && IFDATA(IFcount).FORLOC(2)...
            == IFDATA(OuterIfCount).FORLOC(2)
          % Var is supposed to return as outer overmap - this means that it
          % is not used again within the parent branch and therefore, this
          % branch. In this case, we can just share the same overmap.
          OVERMAP.IF(OVERMAP.IF == OverLoc) = OuterOverLoc;
          IfOvermap(Vcount).Outer = [OuterOverLoc 0];
          IFDATA(OuterIfCount).OVERMAP(OuterOverCount).Inner(end) = [];
          OverCount               = OverCount - 1;
          OverLoc                 = OuterOverLoc;
        elseif ReturnCount == LastOp
          % Is the last variable within the block - see if it gets used
          % again with its branch.
          ReturnFlag = 1;
          if any(IFDATA(IFcount).BROS(BroCount).INNERDEP == LastOp)
            % This variable is used in this branch after its last assignment
            ReturnFlag = 0;
          elseif BroCount ~= nBros
            % This isnt the last branch, see if any of the following
            % branches are dependent upon this variable (i.e. they will
            % change the variable)
            ReturnFlag = 1;
            for Bro2Count = BroCount+1:nBros
              if any(NAMELOCS(IFDATA(IFcount).BROS(Bro2Count).PRIORDEP,1) == NameLoc)
                ReturnFlag = 0;
                break
              end
            end
          end
          if ReturnFlag
            % We can return the overmap after the last union is performed
            RETURN(LastOp,1) = 1;
          else
            % We cant return the overmap after the last union is performed
            % - must return it within the last IfIterEnd
            IfOvermap(Vcount).Return   = NameLoc;
            if RemapFlag
              % We must do the re-mapping within the brother IfIterEnd
              SaveCount                  = SaveCount+1;
              SAVE.IF(LastOp,2)          = SaveCount;
              BROREMAPS{BroCount}(:,Vcount) = [SaveCount; LastOp];
            end
          end
        elseif RemapFlag
          % Is not the last time the var is assigned - see if need to remap
          % within brother IfIterEnd
          if any(IFDATA(IFcount).BROS(BroCount).INNERDEP == LastOp)
            % This variable is used in this branch after its last assignment
            SaveCount                  = SaveCount+1;
            SAVE.IF(LastOp,2)          = SaveCount;
            BROREMAPS{BroCount}(:,Vcount) = [SaveCount; LastOp];
          end
        end
        
      else
        % This variable wasnt assigned within this conditional block, need
        % to use the variable in its state prior to the IF statement.
        % Flag these as "saves" as well
        PrevFlag   = 1;
        PrevCounts = 1:Start-1;
        LastPrevOp = PrevCounts(find(NAMELOCS(PrevCounts,1) == NameLoc,1,'last'));
        if ~isempty(LastPrevOp)
          SaveLoc = SAVE.IF(LastPrevOp,1);
          if ~SaveLoc
            SaveCount = SaveCount+1;
            SaveLoc   = SaveCount;
            SAVE.IF(LastPrevOp,1) = SaveLoc;
          end
          IfOvermap(Vcount).Prior = LastPrevOp;
          OVERMAP.IF(LastPrevOp,4) = OverLoc;
          
          % Check to see if need to do remap within IfIterEnd
          if RemapFlag
            BROREMAPS{BroCount}(:,Vcount) = [Vcount; 0];
          end
        else
          % Wasnt defined prior to if statement either
          PrevFlag = 0;
        end
      end
    end
    
    if ~IFDATA(IFcount).ELSEFLAG && ~PrevFlag
      % Need to check for any variables coming into the IF statement as
      % well
      PrevFlag   = 1;
      PrevCounts = 1:Start-1;
      LastPrevOp = PrevCounts(find(NAMELOCS(PrevCounts,1) == NameLoc,1,'last'));
      OVERMAP.IF(LastPrevOp,4) = OverLoc;
      if ~isempty(LastPrevOp)
        SaveLoc = SAVE.IF(LastPrevOp,1);
        if ~SaveLoc
          SaveCount = SaveCount+1;
          SaveLoc = SaveCount;
          SAVE.IF(LastPrevOp,1) = SaveLoc;
        end
        IfOvermap(Vcount).Prior = LastPrevOp;
      end
    end

    if PrevFlag && INNERLOC && (FunID == 1 || INNERLOC > 1)
      % In a FOR loop - we need to see if the variable coming into the IF
      % statement is different on iteration 1 than on subsequent
      % iterations of a FOR loop, need to do this for all loops which we
      % are embedded within
      for ForLoc = ForLoops
        LoopCounts = FORDATA(ForLoc).START:FORDATA(ForLoc).END;
        LastLoopOp = LoopCounts(find(NAMELOCS(LoopCounts,1) == NameLoc,1,'last'));
        if ~isempty(LastLoopOp)
          SAVE.IF(LastLoopOp,1) = SAVE.IF(LastPrevOp,1);
          if UNROLL && ~any(PRIORDEP == IfOvermap(Vcount).Prior)
            PRIORDEP(end+1) = IfOvermap(Vcount).Prior; %#ok<AGROW>
          end
        end
      end
    end
    
    % -------------- Check for Conditional inside of Loop --------------- %
    if OUTERLOC && ~UNROLL
      % This cond block is embedded within a loop - need to edit the loop
      % overmapping scheme such that all objects belonging to a conditional
      % overmap share the same loop overmap
      LoopStart  = FORDATA(OUTERLOC).START;
      LoopEnd    = FORDATA(OUTERLOC).END;
      LoopCounts = LoopStart:LoopEnd;
      % Need to figure out which loop overmap spot to use
      if OuterOverLoc
        % This cond overmap belongs to a parent cond overmap - use the
        % LastIfCount as this obj will also belong to the parent cond
        % overmap
        ForOverLoc = OVERMAP.FOR(LastIfCount,1);
      elseif PrevFlag && LastPrevOp < LoopStart && OVERMAP.FOR(LastPrevOp,2)
        % obj prior to cond block belongs to cond overmap and is defined
        % prior to loop and belongs to loop overmap, use this loc
        ForOverLoc = OVERMAP.FOR(LastPrevOp,2);
      elseif PrevFlag && LastPrevOp >= LoopStart
        % obj prior to cond block belongs to cond overmap and is defined
        % within loop, use this loc
        ForOverLoc = OVERMAP.FOR(LastPrevOp,1);
      else
        % just use the last object belonging to this cond overmap
        ForOverLoc = OVERMAP.FOR(LastIfCount,1);
      end
      IfOverCounts = AllVarCounts(OVERMAP.IF(AllVarCounts,1)==OverLoc);
      for IOcount = IfOverCounts
        OldForOverLoc = OVERMAP.FOR(IOcount,1);
        OVERMAP.FOR(OVERMAP.FOR == OldForOverLoc) = ForOverLoc;
      end
      
      if PrevFlag
        if LastPrevOp < LoopStart
          % Var coming into this cond block and shares this overmap is
          % defined prior to the outer loop
          if ~OVERMAP.FOR(LastPrevOp,2)
            % Is not marked as a loop overmap prior to loop - do this
            OVERMAP.FOR(LastPrevOp,2) = ForOverLoc;
            SAVE.FOR(LastPrevOp,1)    = SAVE.IF(LastPrevOp,1);
            FORDATA(OUTERLOC).PREVOVERMAP(end+1) = LastPrevOp;
          end
        end
        % Need to make all vars which can be inputs to the block share
        % the same overmap - have already searched for these and they
        % will be marked by the SAVE.IF(i,1)
        IfOverCounts = LoopCounts(SAVE.IF(LoopCounts,1) == SAVE.IF(LastPrevOp,1));
        for IOcount = IfOverCounts
          OldForOverLoc = OVERMAP.FOR(IOcount,1);
          OVERMAP.FOR(OVERMAP.FOR == OldForOverLoc) = ForOverLoc;
        end
      end
    end
  end

  % Mark all variables which this block is dependent upon (but get changed
  % within the block) to be saved.
  PRIORDEP = unique(PRIORDEP);
  if ~isempty(PRIORDEP)
    for VarCount = PRIORDEP
      SaveLoc = SAVE.IF(VarCount,1);
      if ~SaveLoc
        SaveCount = SaveCount+1;
        SaveLoc   = SaveCount;
        SAVE.IF(VarCount,1) = SaveLoc;
      end
      if INNERLOC && (FunID == 1 || INNERLOC > 1)
        % In a FOR loop - we need to see if the variable coming into the IF
        % statement is different on iteration 1 than on subsequent
        % iterations of a FOR loop, need to do this for all loops which we
        % are embedded within
        NameLoc = NAMELOCS(VarCount,1);
        for ForLoc = ForLoops
          LoopCounts = FORDATA(ForLoc).START:FORDATA(ForLoc).END;
          LastLoopOp = LoopCounts(find(NAMELOCS(LoopCounts,1) == NameLoc,1,'last'));
          SaveLoc2 = SAVE.IF(LastLoopOp,1);
          if SaveLoc2 && SaveLoc2 ~= SaveLoc
            error('not sure that this can happen')
          end
          SAVE.IF(LastLoopOp,1) = SaveLoc;
        end
      end
    end
  end

  % Re-Maps that need to be made at the end of each branch
  for BroCount = 1:nBros
    if RemapFlag
      IFDATA(IFcount).BROS(BroCount).REMAPS = ...
        BROREMAPS{BroCount}(:,logical( BROREMAPS{BroCount}(1,:)) );
    else
      IFDATA(IFcount).BROS(BroCount).REMAPS = [];
    end
  end
  % Overmapping Scheme to let routines know what needs to be changed prior
  % to and at the end of the block
  IFDATA(IFcount).OVERMAP = IfOvermap;
  IFDATA(IFcount).PRIORDEP.Locs = PRIORDEP;
end
if ~isempty(OVERMAP.IF)
  IFOVERLOCS = OVERMAP.IF(:,5);
  OVERMAP.IF = OVERMAP.IF(:,[1 4]);
end

%% ~~~~~~~~~~~~~~~~~~~~~~ BREAKS/CONTINUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
if ~UNROLL
  for ContCount = 1:size(CONTLOCS,1)+size(BREAKLOCS,1)
    
    if ContCount <= size(CONTLOCS,1)
      IfLoc   = CONTLOCS(ContCount,1);
      BroLoc  = CONTLOCS(ContCount,2);
      ForLoc  = CONTLOCS(ContCount,3);
      oForLoc = CONTLOCS(ContCount,4);
    else
      IfLoc   = BREAKLOCS(ContCount,1);
      BroLoc  = BREAKLOCS(ContCount,2);
      ForLoc  = BREAKLOCS(ContCount,3);
      oForLoc = BREAKLOCS(ContCount,4);
    end
    
    
    bStart = IFDATA(IfLoc).BROS(BroLoc).START;
    bEnd   = IFDATA(IfLoc).BROS(BroLoc).END;
    iStart = IFDATA(IfLoc).BROS(1).START;
    fStart  = FORDATA(ForLoc).START;   % START OF INNERMOST FOR LOOP
    fEnd    = FORDATA(ForLoc).END;     % END OF INNERMOST FOR LOOP
    ofStart = FORDATA(oForLoc).START;  % START OF OUTERMOST FOR LOOP
    
    cVarCounts = [fStart:iStart-1 bStart:bEnd];
    % active variables at time of break/continue firing
    pVarCounts = 1:fStart-1;
    % active variables prior to the loop
    if ContCount <= size(CONTLOCS,1)
      % CONTINUE statement overmapping - for each variable active at time of
      % the statement:
      % IF
      % 1. the variable gets a different assignment if the continue statement
      %    doesnt fire
      % 2. the variable occurs prior to the loop
      % 3. the loop is dependent upon the variable outside of the loop
      % THEN the active instance at the time of the continue statement firing
      % must share a loop overmap with the instances from 1. and 2.
      %
      % ALSO, after each loop iteration, we must return the UNION of the
      % continue statement firing and it not firing.
      
      % So, the CONTINUE OVERMAP gets cleared at the start of each iteration,
      % gets built by the union of the active workspace at the time of
      % continue firing and the active workspace at the end of the loop, on
      % iterations 1:N-1 where N is the number of iterations in the loop. On
      % the last iteration, the CONTINUE statement is treated as a BREAK
      % statement - see below.
      cAsgnCounts = cVarCounts(logical(NAMELOCS(cVarCounts,1)));
      cNameLocs   = unique(NAMELOCS(cAsgnCounts,1));
      % if we need to make this connection, then there will have been
      % something assigned prior to the loop to make a similar connection
      for Ncount = 1:length(cNameLocs)
        NameLoc = cNameLocs(Ncount);
        cVarLoc = cVarCounts(find(NAMELOCS(cVarCounts,1)==NameLoc,1,'last'));
        pVarLoc = pVarCounts(find(NAMELOCS(pVarCounts,1)==NameLoc,1,'last'));
        
        if ~isempty(pVarLoc) && LASTOCC(pVarLoc) >= fStart && LASTOCC(pVarLoc) <= fEnd
          % 2. and 3. are satisfied
          if cVarLoc >= bStart
            % Variable is defined within conditional branch that contains
            % continue statement - variables prior to the conditional set
            % belong to the loop
            fVarCounts = [fStart:bStart-1 bEnd+1:fEnd];
          else
            % Variable is NOT defined within the conditional branch that
            % contains continue statement - variables prior to the conditional
            % set DO NOT belong to the loop
            fVarCounts = [iStart:bStart-1 bEnd+1:fEnd];
          end
          fVarLoc = fVarCounts(find(NAMELOCS(fVarCounts,1)==NameLoc,1,'last'));
          if ~isempty(fVarLoc)
            % 1. is satisfied
            if pVarLoc < ofStart
              pOverLoc = OVERMAP.FOR(pVarLoc,2);
            else
              pOverLoc = OVERMAP.FOR(pVarLoc,1);
            end
            cOverLoc = OVERMAP.FOR(cVarLoc,1);
            fOverLoc = OVERMAP.FOR(fVarLoc,1);
            % Set all loop overmaps to share the same vars
            OVERMAP.FOR(OVERMAP.FOR == cOverLoc) = pOverLoc;
            OVERMAP.FOR(OVERMAP.FOR == fOverLoc) = pOverLoc;
            
            % Set up CONTINUE overmap (to be enforced after each iteration
            % 1:N-1 (during the Overmapping Evaluation)
            fOverLoc = OVERMAP.CONT(fVarLoc);
            cOverLoc = OVERMAP.CONT(cVarLoc);
            if fOverLoc && cOverLoc
              OVERMAP.CONT(OVERMAP.CONT == cOverLoc) = fOverLoc;
            elseif fOverLoc
              OVERMAP.CONT(cVarLoc) = fOverLoc;
            elseif cOverLoc
              OVERMAP.CONT(fVarLoc) = cOverLoc;
            else
              OverCount = OverCount+1;
              OVERMAP.CONT([cVarLoc fVarLoc]) = OverCount;
            end
          end
        end
      end
    end
    
    % BREAK statement overmapping - for each variable active at the time of
    % the statement:
    % IF
    % 1. the variable gets a different assignment if the break doesnt fire
    % 2. the variable is used after the loop (either in a parent loop or
    % outside of the loop)
    % THEN the active instance at the time of the break statement firing
    % must share a loop overmap with the instance from 1.
    %
    % ALSO, after all iterations of the loop have been run, we must UNION
    % the variables from each iteration of the BREAK with the result of the
    % last iteration of the loop.
    
    % We also note that a CONTINUE statement on the last iteration of a
    % loop is the same thing as a BREAK statement
    
    % Need to find which loops the INNER loop is embedded within
    ForLoops = zeros(1,ForLoc-oForLoc);
    forLoci  = ForLoc;
    for iForCount = 1:ForLoc-oForLoc
      forLoci = FORDATA(forLoci).PARENTLOC;
      if ~isempty(forLoci)
        ForLoops(iForCount) = forLoci;
      else
        ForLoops = ForLoops(1:iForCount-1);
        break
      end
    end
    if FunID > 1 && FunAsLoopFlag; ForLoops(end) = []; end
    
    cAsgnCounts = cVarCounts(logical(NAMELOCS(cVarCounts,1)));
    cNameLocs   = unique(NAMELOCS(cAsgnCounts,1));
    for Ncount = 1:length(cNameLocs)
      NameLoc = cNameLocs(Ncount);
      cVarLoc = cVarCounts(find(NAMELOCS(cVarCounts,1)==NameLoc,1,'last'));
      if cVarLoc >= bStart
        % Variable is defined within conditional branch that contains
        % break statement - variables prior to the conditional set
        % belong to the loop
        fVarCounts = [fStart:bStart-1 bEnd+1:fEnd];
      else
        % Variable is NOT defined within the conditional branch that
        % contains continue statement - variables prior to the conditional
        % set DO NOT belong to the loop
        fVarCounts = [iStart:bStart-1 bEnd+1:fEnd];
      end
      fVarLoc = fVarCounts(find(NAMELOCS(fVarCounts,1)==NameLoc,1,'last'));
      if ~isempty(fVarLoc)
        % 1. is satisfied
        if cVarLoc > fVarLoc; LastLoc = cVarLoc; else LastLoc = fVarLoc; end
        UnionFlag = 0;
        if LASTOCC(LastLoc) > fEnd
          % the variable gets used after the loop, so we need to union the
          % two.
          UnionFlag = 1;
        else
          % need to check and see if this is a nested loop and the variable
          % gets used in the nested loop.
          for forLoci = 1:ForLoops
            forStarti = FORDATA(forLoci).START;
            forEndi   = FORDATA(forLoci).END;
            if ~any(NAMELOCS(fEnd+1:forEndi,1) == NameLoc)
              % isnt reassigned after the innermost loop
              prevCounts = 1:forStarti-1;
              prevCount  = prevCounts(find(NAMELOCS(prevCounts,1)==NameLoc,1,'last'));
              if ~isempty(prevCount) && LASTOCC(prevCount) > forStarti
                % variable is defined outside of a parent loop, used within
                % the parent loop
                UnionFlag = 1;
                break
              end
            end
          end
        end
        if UnionFlag
          % Make all these variables share the same LOOP OVERMAP
          cOverLoc = OVERMAP.FOR(cVarLoc,1);
          fOverLoc = OVERMAP.FOR(fVarLoc,1);
          if cVarLoc > fVarLoc
            OVERMAP.FOR(OVERMAP.FOR == fOverLoc) = cOverLoc;
          else
            OVERMAP.FOR(OVERMAP.FOR == cOverLoc) = fOverLoc;
          end
          
          % Assign the BREAK OVERMAP - we assign a minus to variables that
          % only get used on the last iteration and a plus to variables that
          % get used on each iteration i.e. BREAK variables get a plus,
          % CONTINUE and neither get a minus
          if ContCount <= size(CONTLOCS,1)
            % Since we do all continues first, we dont need to worry about
            % overriding any plus signs.
            cOverLoc = OVERMAP.BREAK(cVarLoc,1);
            fOverLoc = OVERMAP.BREAK(fVarLoc,1);
            if fOverLoc && cOverLoc
              OVERMAP.BREAK(OVERMAP.BREAK == cOverLoc) = fOverLoc;
            elseif fOverLoc
              OVERMAP.BREAK(cVarLoc) = fOverLoc;
            elseif cOverLoc
              OVERMAP.BREAK(fVarLoc) = cOverLoc;
            else
              OverCount = OverCount+1;
              OVERMAP.BREAK([cVarLoc fVarLoc]) = -OverCount;
            end
          else
            % This is for a BREAK, plus overrides minus, but minus does not
            % override plus
            cOverLoc = OVERMAP.BREAK(cVarLoc,1);
            fOverLoc = OVERMAP.BREAK(fVarLoc,1);
            if fOverLoc && cOverLoc
              if fOverLoc > 0; OverLoc = fOverLoc; else OverLoc = -fOverLoc; end
              % Using fOver, so dont need to look to change any of those
              if cOverLoc > 0
                OVERMAP.BREAK(OVERMAP.BREAK == cOverLoc) = OverLoc;
                OVERMAP.BREAK(OVERMAP.BREAK == -cOverLoc) = -OverLoc;
              else
                OVERMAP.BREAK(OVERMAP.BREAK == -cOverLoc) = OverLoc;
                OVERMAP.BREAK(OVERMAP.BREAK == cOverLoc) = -OverLoc;
                OVERMAP.BREAK(cVarLoc) = OverLoc;
              end
            elseif fOverLoc
              if fOverLoc > 0
                OVERMAP.BREAK(cVarLoc) = fOverLoc;
              else
                OVERMAP.BREAK(cVarLoc) = -fOverLoc;
              end
            elseif cOverLoc
              if cOverLoc > 0
                OVERMAP.BREAK(fVarLoc) = -cOverLoc;
              else
                OVERMAP.BREAK(cVarLoc) = -cOverLoc;
                OVERMAP.BREAK(fVarLoc) = cOverLoc;
              end
            else
              OverCount = OverCount+1;
              OVERMAP.BREAK(cVarLoc) = OverCount;
              OVERMAP.BREAK(fVarLoc) = -OverCount;
            end
          end
        end
      end
    end
  end
end
if ~isempty(OVERMAP.IF)
  OVERMAP.IF = [OVERMAP.IF(:,1) IFOVERLOCS];
end
VARSTORAGE.OVERMAP = cell(OverCount,1);
VARSTORAGE.SAVE    = cell(SaveCount,1);
VARINFO.OVERMAP = OVERMAP;
VARINFO.SAVE    = SAVE;
VARINFO.RETURN  = RETURN;
end