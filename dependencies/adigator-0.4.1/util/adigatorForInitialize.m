function [myLoopVar, outEvalStr, outEvalVar] = adigatorForInitialize(ForCount,UserLoopVar)
% This transformation routine is called prior to the evaluation of any FOR
% loop in the intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
if ADIGATOR.OPTIONS.PREALLOCATE
  % --------------------------------------------------------------------- %
  %                    Cell/Struct Pre-Allocation Run                     %
  % --------------------------------------------------------------------- %
  myLoopVar = 1:size(UserLoopVar,2);
  outEvalStr = [];
  outEvalVar = [];
elseif ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  ADIGATORFORDATA(ForCount).COUNTNAME = sprintf('cadaforcount%1.0d',ForCount);
  
  ADIGATORFORDATA(ForCount).START     = ADIGATOR.VARINFO.COUNT;
  ADIGATORFORDATA(ForCount).PARENTLOC = ADIGATOR.FORINFO.INNERLOC;
  
  if ADIGATOR.FORINFO.FLAG
    ADIGATORFORDATA(ADIGATOR.FORINFO.INNERLOC).CHILDLOCS(end+1,1) = ForCount;
  elseif ~ADIGATOR.FORINFO.OUTERLOC
    ADIGATOR.FORINFO.OUTERLOC = ForCount;
  end
  ADIGATOR.FORINFO.INNERLOC = ForCount;
  ADIGATORFORDATA(ForCount).CHILDLOCS = zeros(0,1);
  if ~ADIGATOR.OPTIONS.UNROLL
    ADIGATOR.FORINFO.FLAG = 1;
  end
  outEvalStr = [];
  outEvalVar = [];
  myLoopVar = 0;
  if ADIGATOR.OPTIONS.PREALLOCATE
    % We need to pre-allocate...
    if isa(UserLoopVar,'cada')
      ForLength = UserLoopVar.func.size(2);
      func.size = [1 1];
      func.loopvar = ForLength;
      myLoopVar = cada(0,func,[]);
    end
  end
  
elseif ADIGATOR.OPTIONS.UNROLL
  % --------------------------------------------------------------------- %
  %            Unrolling Loop - Overmapping or Printing Eval              %
  % --------------------------------------------------------------------- %
  if ADIGATOR.EMPTYFLAG
    % Dont run this loop
    myLoopVar = [];
  elseif isa(UserLoopVar,'cada')
    ForLength = UserLoopVar.func.size(2);
    if isinf(ForLength)
      error('Cannot loop over vectorized dimension')
    end
    myLoopVar = 1:ForLength;
  else
    myLoopVar = UserLoopVar;
    ForLength = length(myLoopVar);
  end
  % ---------------------- Set For Lengths ---------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;

  if ParentLoc
    ParentIter = ADIGATORFORDATA(ParentLoc).COUNT.ITERATION;
    if ParentIter == 1
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
    else
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(ParentIter) = ForLength;
    end
  else
    ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
  end
  
  if isempty(myLoopVar)
    % Dont need to run loop
    ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).END+1;
    ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  else
    ADIGATOR.FORINFO.INNERLOC = ForCount;
    if ~ParentLoc
      ADIGATOR.FORINFO.OUTERLOC = ForCount;
    end
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if ADIGATOR.EMPTYFLAG
    % Dont run this loop
    ForLength = 0;
    myLoopVar = [];
  elseif isa(UserLoopVar,'cada')
    ForLength = UserLoopVar.func.size(2);
    if isinf(ForLength)
      error('Cannot loop over vectorized dimension')
    end
    myLoopVar = 1:ForLength;
  else
    myLoopVar = UserLoopVar;
    ForLength = length(myLoopVar);
  end
  % ---------------------- Set For Lengths ---------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;
  if ParentLoc
    ParentIter = ADIGATORFORDATA(ParentLoc).COUNT.ITERATION;
    if ParentIter == 1
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
    else
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(ParentIter) = ForLength;
    end
  else
    ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
  end
  
  if isempty(myLoopVar)
    % Dont ned to run this loop. Set exit sequence
    ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).END+1;
    ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  else
    % -------------- Overmap Variables Coming into Loop ----------------- %
    ADIGATOR.FORINFO.INNERLOC = ForCount;
    if ~ParentLoc
      ADIGATOR.FORINFO.OUTERLOC = ForCount;
    end
    for VarCount = ADIGATORFORDATA(ForCount).PREVOVERMAP
      SaveLoc = ADIGATOR.VARINFO.SAVE.FOR(VarCount,1);
      OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,2);
      PrevVar = ADIGATORVARIABLESTORAGE.SAVE{SaveLoc};
      OverVar = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
      if isa(OverVar,'cada')
        OverVar = cadaUnionVars(PrevVar,OverVar);
      else
        OverVar = PrevVar;
      end
      ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc} = OverVar;
    end
    ADIGATOR.FORINFO.FLAG = 1;
  end
  outEvalStr = [];
  outEvalVar = [];
  
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;
  ADIGATOR.FORINFO.INNERLOC = ForCount;
  ADIGATOR.FORINFO.FLAG     = 1;
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  if ~ParentLoc
    % Is an Outer Loop
    ADIGATOR.FORINFO.OUTERLOC = ForCount;
    % -------------- Remap Variables Coming Into Loop ------------------- %
    nOutEval = length(ADIGATORFORDATA(ForCount).PREVOVERMAP);
    outEvalStr = cell(nOutEval,1);
    outEvalVar = cell(nOutEval,1);
    for Vcount = 1:nOutEval
      % Do Re-Map
      VarCount = ADIGATORFORDATA(ForCount).PREVOVERMAP(Vcount);
      NameLoc  = ADIGATOR.VARINFO.NAMELOCS(VarCount,1);
      VarStr   = ADIGATOR.VARINFO.NAMES{NameLoc};
      OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,2);
      PrevVar  = ADIGATORVARIABLESTORAGE.SAVE{ADIGATOR.VARINFO.SAVE.FOR(VarCount,1)};
      if OverLoc
        OverVar  = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
        outEvalVar{Vcount} = cadaPrintReMap(PrevVar,OverVar,VarCount);
      else
        outEvalVar{Vcount} = PrevVar;
      end
      outEvalStr{Vcount} = sprintf([VarStr,' = adigatorForEvalVar{%1.0d};'],Vcount);
      % Check for case when conditional set needs to return a saved
      % variable.
      if ~isempty(ADIGATOR.VARINFO.SAVE.IF)
        IfSaveLoc = ADIGATOR.VARINFO.SAVE.IF(VarCount,1);
        if IfSaveLoc
          ADIGATORVARIABLESTORAGE.SAVE{IfSaveLoc} = outEvalVar{Vcount};
        end
      end
    end
    
    % ---------------------- Print the Loop ----------------------------- %
    if ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
        ' = 1:%1.0d\n'],ADIGATORFORDATA(ForCount).MAXLENGTH);
    end
  else
    % Is an Inner Loop
    % -------Print Out any Organizational Function Data Reshape/Refs----- %
    OrgFuncs = {'SUBSREF','SUBSASGN','SPARSE','NONZEROS','HORZCAT',...
      'VERTCAT','TRANSPOSE','REPMAT','RESHAPE','SIZE'};
    for OFcount = 1:10
      DataStructure = ADIGATORFORDATA(ForCount).(OrgFuncs{OFcount});
      % do .INDICES fields first
      if strcmp(OrgFuncs{OFcount},'HORZCAT') ||...
          strcmp(OrgFuncs{OFcount},'VERTCAT')
        % horzcat and vertcat have added dimension
        for OF2count = 1:length(DataStructure)
          StringArray = DataStructure(OF2count).INDICES;
          for Icount = 1:size(StringArray,1)
            for Jcount = 1:size(StringArray,2)
              LHSstr = StringArray{Icount,Jcount,1};
              RHSstr = StringArray{Icount,Jcount,2};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      elseif isfield(DataStructure,'INDICES')
        % will get the other ones
        for OF2count = 1:length(DataStructure)
          StringArray  = DataStructure(OF2count).INDICES;
          SubsAsgnFlag = size(StringArray,2)>3;
          for Icount = 1:size(StringArray,1)
            LHSstr = StringArray{Icount,1};
            RHSstr = StringArray{Icount,2};
            if ~isempty(LHSstr) && ~isempty(RHSstr)
              fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
            end
            if SubsAsgnFlag % checks for subsasgn, he has extra
              LHSstr = StringArray{Icount,4};
              RHSstr = StringArray{Icount,5};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      end
      if isfield(DataStructure,'SIZES')
        % do .SIZES field next - all are straightforward
        for OF2count = 1:length(DataStructure)
          StringArray = DataStructure(OF2count).SIZES;
          if iscell(StringArray)
            for Icount = 1:size(StringArray,1)
              LHSstr = StringArray{Icount,1};
              RHSstr = StringArray{Icount,2};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      end
    end

    % print out any FOR loop changing size stuff
    for Fcount = 1:length(ADIGATORFORDATA(ForCount).FOR)
      LHSstr = ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS{1};
      RHSstr = ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS{2};
      if ~isempty(LHSstr) && ~isempty(RHSstr)
        fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
      end
    end
    
    % -------------------- Print the Loop ------------------------------- %
    if ~ADIGATOR.EMPTYFLAG
      if ADIGATOR.DERNUMBER == 1
        LoopVarStr = ADIGATORFORDATA(ForCount).FOR(1).LENGTHS{1};
        if ~isempty(LoopVarStr)
          % Size of loop is dependent upon an outer loop
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:',LoopVarStr,'\n']);
        else
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:%1.0d\n'],ADIGATORFORDATA(ForCount).MAXLENGTH);
        end
      else
        if isa(UserLoopVar,'cada')
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:',UserLoopVar.func.name,'\n']);
        else
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:%1.0d\n'],length(UserLoopVar));
        end
      end
    end
    outEvalStr = [];
    outEvalVar = [];
  end
  ADIGATOR.PRINT.INDENT = [indent,'    '];
  myLoopVar = 1;
end

if ADIGATOR.OPTIONS.UNROLL
  ADIGATOR.FORINFO.FLAG = 0;
end