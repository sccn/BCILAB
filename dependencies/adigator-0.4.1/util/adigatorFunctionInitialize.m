function [flag, FunctionInfo, Outputs, varargout] = adigatorFunctionInitialize(FunID,FunctionInfo,Inputs)
% function [flag, Outputs] = adigatorFunctionInitialize(FunID,Inputs)
% This function is called from the evaluating files when the main function
% or any sub-function is initially called.
% ------------------------ Input Information ---------------------------- %
% FunID :  function ID - what element of FunctionInfo this func is
% Inputs:  cell array of the inputs to the function
% ----------------------- Output Information ---------------------------- %
% flag:     this is set to 0 or 1, it is set to 0 if the function is to be
%           evaluated, 1 if the function is not to be evaluated.
% Outputs:  these are either the actual inputs to the function or the
%           outputs of the function. If flag = 1, then these are the
%           outputs (the function didnt need to be run), if flag = 0, then
%           these are the inputs.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR  ADIGATORDATA ADIGATORFORDATA ADIGATORVARIABLESTORAGE

NUMinputs = length(Inputs);
if any(ADIGATOR.FILE.PARENTID == FunID)
  error('Functions cannot call themselves from within their methods.');
end
% ----------------------------------------------------------------------- %
%                            PARSE INPUTS                                 %
% ----------------------------------------------------------------------- %
PDflag = FunctionInfo(FunID).DERNUMBER > 1 && ADIGATOR.DERNUMBER == 1;
if (ADIGATOR.RUNFLAG == 2 && ~FunctionInfo(FunID).Iteration.CallCount) ||...
    (ADIGATOR.RUNFLAG && ADIGATOR.OPTIONS.UNROLL)
  % Printing this function, set these so we can get proper naming
  if ADIGATOR.OPTIONS.UNROLL && ~isempty(ADIGATOR.FILE.PARENTID)
    FunctionInfo(ADIGATOR.FILE.PARENTID(end)).VARINFO = ADIGATOR.VARINFO;
  end
  ADIGATOR.VARINFO         = FunctionInfo(FunID).VARINFO;
  ADIGATOR.DERNUMBER       = FunctionInfo(FunID).DERNUMBER;
  ADIGATOR.FILE.FUNID      = FunID;
end

if ADIGATOR.OPTIONS.PREALLOCATE
  % Running this numerically to preallocate cells/structures
  ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;
  Outputs = Inputs; flag = 0;
  if FunctionInfo(FunID).DERNUMBER > 1
    varargout{1} = FunctionInfo(FunID).PreviousDerivData;
  else
    varargout = cell(1);
  end
  return
end

NumInVars   = 0;
InVars      = cell(NUMinputs,1);
InVarStrs   = cell(NUMinputs,1);
InOps       = zeros(NUMinputs,1);
for Icount = 1:NUMinputs
  %% ~~~~~~~~~~~~~~~~~~~~~~~~~ PARSE INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  % We are looking for:
  %   1. Overloaded Objects - do nothing with these
  %   2. Numeric Objects - turn these into overloaded objects (according
  %     to what option the user has set, default is that they hold the
  %     same sparsity pattern.
  %   3. Structures/Cells - for these we need to scan them until we get
  %     to either an Overloaded Object or a Numeric Object.
  CurVar    = Inputs{Icount};
  CurVarStr = FunctionInfo(FunID).Input.Names{Icount};
  [Inputs{Icount},InVars,InVarStrs,InOps,NumInVars] =...
    FindInputVars(CurVar,CurVarStr,InVars,InVarStrs,InOps,NumInVars,PDflag);
end
InVars    = InVars(1:NumInVars);
InVarStrs = InVarStrs(1:NumInVars);
InOps     = InOps(1:NumInVars);
% InVars are now all of the overloaded variables which the function starts
% with. InVarStrs are their string names (with structure and cell reference
% extensions).

if PDflag && ~ADIGATOR.RUNFLAG
  ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC >= min(InOps),2) = -FunID;
  PrevAux = FunctionInfo(FunID).PreviousDerivData.AuxInputs;
  AuxInputFlags = zeros(size(PrevAux));
  for PAcount = 1:length(PrevAux)
    InAuxLoc = strcmp(PrevAux{PAcount},InVarStrs);
    InAuxOp  = InOps(InAuxLoc);
    if ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC == InAuxOp,3) == -Inf
      % These are also Aux inputs to the main function
      ADIGATOR.VARINFO.NAMELOCS(InAuxOp,3) = -Inf;
      AuxInputFlags(PAcount) = 1;
    end
    % remove the -FunID in column 2 for aux inputs
    ADIGATOR.VARINFO.NAMELOCS(InAuxOp,2) = 0;
    ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC == InAuxOp,2) = 0;
  end
  FunctionInfo(FunID).PreviousDerivData.AuxInputFlags = AuxInputFlags;
end
if FunctionInfo(FunID).DERNUMBER > 1 && ...
    strcmp(FunctionInfo(FunID).PreviousDerivData.Derivative...
    (FunctionInfo(FunID).DERNUMBER-1).FunType,'Main')
  %% ~~~~~~~~~~~~~~~~ CHECK INPUTS IF 2ND or HIGHER DERIV ~~~~~~~~~~~~~~ %%
  if ~ADIGATOR.RUNFLAG
    % Make naming scheme of inputs as if calling function was 2nd deriv
  else
    InputChecks = FunctionInfo(FunID).PreviousDerivData...
      .Derivative(FunctionInfo(FunID).DERNUMBER-1).InputChecks;
    for Icount = 1:length(InputChecks)
      derivChecks = InputChecks(Icount).deriv;
      for Dcount = 1:length(derivChecks)
        dername = derivChecks(Dcount).vodname;
        for D2count = 1:ADIGATOR.NVAROFDIFF
          if strcmp(ADIGATOR.VAROFDIFF(D2count).name,dername)
            % Taking another derivative wrt this variable - check that nzlocs
            % are defined as they were previously.
            inputloc = strcmp(InputChecks(Icount).name,InVarStrs);
            oldnzlocs = derivChecks(Dcount).nzlocs;
            newnzlocs = InVars{inputloc}.deriv(D2count).nzlocs;
            if ~isequal(oldnzlocs,newnzlocs)
              error(['input ''',InputChecks(Icount).name,''' should have non-zero ',...
                'derivatives wrt ''',dername,''' in locations defined by ',mat2str(oldnzlocs)]);
            end
          end
        end
      end
    end
  end
end


if ~ADIGATOR.RUNFLAG
  %% ~~~~~~~~~~~~~~~~~~~~~~~~~ EMPTY RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  FunctionInfo = storeOldGlobalParams(FunctionInfo);
  ADIGATOR.FILE.FUNID    = FunID;
  ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
  FunAsLoopFlag = ...
    ((ADIGATOR.FORINFO.FLAG && ~FunctionInfo(FunID).Iteration.CallCount) ||...
      (FunctionInfo(FunID).Iteration.CallCount && ...
      ~FunctionInfo(FunID).FunAsLoopFlag)) && ~ADIGATOR.OPTIONS.UNROLL;
  if FunAsLoopFlag || ~FunctionInfo(FunID).Iteration.CallCount
    % ------------------------ Run Empty -------------------------------- %
    
    FunctionInfo(FunID).FunAsLoopFlag    = FunAsLoopFlag;
    FunctionInfo(FunID).Input.StrucNames = InVarStrs;
    FunctionInfo(FunID).Input.StrucVars  = InVars;
    
    % Set-Up Global Variables
    ADIGATORFORDATA  = FunctionInfo(FunID).FORDATA;
    ADIGATOR.IFDATA  = FunctionInfo(FunID).IFDATA;
    ADIGATOR.VARINFO.NAMELOCS      = zeros(NumInVars,3);
    ADIGATOR.VARINFO.NAMELOCS(:,1) = 1:NumInVars;

    ADIGATOR.VARINFO.NAMELOCS(:,3) = 1;
    ADIGATOR.VARINFO.LASTOCC       = (1:NumInVars).';
    ADIGATOR.VARINFO.NAMES         = InVarStrs;
    ADIGATOR.VARINFO.COUNT         = NumInVars+1;
    
    ADIGATOR.DERNUMBER     = FunctionInfo(FunID).DERNUMBER;
    ADIGATOR.PREOPCOUNT    = NumInVars+1;
    ADIGATOR.EMPTYFLAG     = 1;
    
    ADIGATOR.BREAKLOCS = zeros(0,4);
    ADIGATOR.CONTLOCS  = zeros(0,4);
    ADIGATOR.ERRORLOCS = zeros(0,2);
    ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;
    
    % Check for numeric inputs and set them to be named differently.
    if FunID == 1
      for Icount = 1:NumInVars
        derflag = cadaCheckForDerivs(InVars{Icount});
        if ~derflag
          ADIGATOR.VARINFO.NAMELOCS(Icount,3) = -Inf;
        end
      end
    else
      ParentID = ADIGATOR.FILE.PARENTID(end-1);
      ParLASTOCC = FunctionInfo(ParentID).VARINFO.LASTOCC;
      ParNAMELOCS = FunctionInfo(ParentID).VARINFO.NAMELOCS;
      for Icount = 1:NumInVars
        if ParNAMELOCS(ParLASTOCC == InOps(Icount),3) == -Inf
          ADIGATOR.VARINFO.NAMELOCS(Icount,3) = -Inf;
        elseif ParNAMELOCS(ParLASTOCC == InOps(Icount),3) == Inf
          ADIGATOR.VARINFO.NAMELOCS(Icount,3) = Inf;
        end
      end
    end
    
    if FunAsLoopFlag
      % Setting this up as the sub-function being the outermost loop
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 1;
      ADIGATOR.FORINFO.OUTERLOC = 1;  ADIGATOR.FORINFO.FLAG = 1;
      ADIGATOR.FORINFO.INNERLOC = 1;  ADIGATOR.FORINFO.FUNLOOP = 1;
      if FunctionInfo(FunID).Iteration.CallCount
        % This was already ran once (under no loop assumption) - clear the
        % collected loop data so we can do it again
        NumLoops = length(ADIGATORFORDATA);
        ADIGATORFORDATA = struct('START',cell(0,1),'END',[],...
          'COUNTNAME',[],'PREVOVERMAP',[],'COUNT',[],'FOR',[],'SUBSREF',...
          [],'SUBSASGN',[],'SPARSE',[],'NONZEROS',[],'HORZCAT',[],'VERTCAT',...
          [],'TRANSPOSE',[],'REPMAT',[],'RESHAPE',[],'SIZE',[],'OTHER',[],...
          'CONTRESTORE',[],'PARENTLOC',[]);
        ADIGATORFORDATA(NumLoops,1).COUNT = [];
        FunctionInfo(FunID).Iteration.CallCount = 0;
      end
      ADIGATORFORDATA(1).COUNTNAME = sprintf('adigator%1.0dfuncount',...
        ADIGATOR.DERNUMBER);
      ADIGATORFORDATA(1).FOR(1).LENGTHS = 1;
      ADIGATORFORDATA(1).CHILDLOCS = zeros(0,1);
      ADIGATORFORDATA(1).PARENTLOC = 0;
      ADIGATORFORDATA(1).START     = 1;
      adigatorForIterStart(1,0);
    else
      % Not making the sub-function a loop
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 0;
      ADIGATOR.FORINFO.INNERLOC = 0;  ADIGATOR.FORINFO.OUTERLOC = 0;
      ADIGATOR.FORINFO.FLAG     = 0;  ADIGATOR.FORINFO.FUNLOOP  = 0;
    end
    Outputs = Inputs; flag = 0;
  else
    % --------------------- Called Again - Return Empty ----------------- %
    Outputs = FunctionInfo(FunID).Output.EmptyVars;
    [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
    flag = 1;
  end
  
elseif ADIGATOR.OPTIONS.UNROLL && FunID > 1
  %% ~~~~~~~~~~~~~~~~~~~ UNROLLING SUB-FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~ %%
  ParentID = ADIGATOR.FILE.PARENTID(end);
  ADIGATOR.VARINFO   = FunctionInfo(ParentID).VARINFO;
  FunctionInfo = storeOldGlobalParams(FunctionInfo);

  ADIGATOR.FILE.FUNID    = FunID;
  ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];

  if ADIGATOR.EMPTYFLAG
    % ------------------------ Empty Run -------------------------------- %
    Outputs = FunctionInfo(FunID).Output.EmptyVars;
    [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
    flag = 1;
  else
    % ----------------- Determine if Need to Print ---------------------- %
    if ~FunctionInfo(FunID).Iteration.CallCount
      SameFlag = 0;
      FunctionInfo(FunID).Iteration.IterCount = 0;
    else
      if ~isequal(FunctionInfo(FunID).Input.StrucNames,InVarStrs)
        % Check and make sure the input structure is the same
        error(['Currently if using multiple calls to a function, the',...
          ' function must always have the same I/O structure, with the',...
          ' same Structure/Cell patterns of the inputs/outputs.'])
      end
      % All we really care about is the input structure, we are unrolling
      % so could care less about where it gets called from.
      for ITcount = 1:FunctionInfo(FunID).Iteration.IterCount
        SameFlag = ITcount;
        for Icount = 1:NumInVars
          InVar1 = InVars{Icount};
          InVar2 = FunctionInfo(FunID).Input.StrucVars{Icount,ITcount};
          InVar1.func.name =[]; InVar2.func.name = [];
          if ~isequal(InVar1.func,InVar2.func)
            SameFlag = 0; break
          else
            for Vcount = 1:ADIGATOR.NVAROFDIFF
              if ~isequal(InVar1.deriv(Vcount).nzlocs,...
                  InVar2.deriv(Vcount).nzlocs)
                SameFlag = 0; break
              end
            end
          end
        end
        if SameFlag; break; end
      end
    end
    CallCount = FunctionInfo(FunID).Iteration.CallCount + 1;
    FunctionInfo(FunID).Iteration.CallCount = CallCount;
    if ~SameFlag
      IterCount = FunctionInfo(FunID).Iteration.IterCount + 1;
      FunctionInfo(FunID).Iteration.IterCount = IterCount;
      CurrentIter = IterCount;
    else
      CurrentIter = SameFlag;
    end
    if ADIGATOR.RUNFLAG == 2 && FunID > 1
      % -------------- Print Input in Calling Func File ----------------- %
      fid = ADIGATOR.PRINT.FID;
      indent = ADIGATOR.PRINT.INDENT;
      InCallStr = cell(1,NUMinputs);
      for Icount = 1:NUMinputs
        InCallStr{Icount} = sprintf('cadainput%1.0d,',Icount);
      end
      InCallStr = cell2mat(InCallStr);
      InCallStr = ['(',InCallStr(1:end-1),')'];

      NUMoutputs = length(FunctionInfo(FunID).Output.Names);
      if NUMoutputs == 1
        OutCallStr = 'cadaoutput1';
      else
        OutCallStr = cell(1,NUMoutputs);
        for Ocount = 1:NUMoutputs
          OutCallStr{Ocount} = sprintf('cadaoutput%1.0d,',Ocount);
        end
        OutCallStr = cell2mat(OutCallStr);
        OutCallStr = ['[',OutCallStr(1:end-1),']'];
      end
      if FunctionInfo(FunID).DERNUMBER > 1
        FunctionStr = FunctionInfo(FunID).File.Name;
        if length(FunctionStr) > 8 && strcmp(FunctionStr(1:7),'ADiGator')
          FunctionStr = ['ADiGator_',FunctionStr(8:end)];
        end
      else
        FunctionStr = sprintf(['ADiGator_',FunctionInfo(FunID).File.Name,'%1.0f'],CurrentIter);
      end
      fprintf(fid,[indent,OutCallStr,' = ',FunctionStr,InCallStr,';\n']);
      if ADIGATOR.OPTIONS.COMMENTS
        fprintf(fid,[indent,'%% Call to function: ',FunctionInfo(FunID).File.Name,'\n']);
      end
    end

    if ~SameFlag
      % --------------------- Run/Print This ---------------------------- %
      FunctionInfo(FunID).Iteration.IterID(CallCount)  = CurrentIter;
      FunctionInfo(FunID).Input.StrucVars(:,IterCount) = InVars;
      
      ADIGATOR.DERNUMBER     = FunctionInfo(FunID).DERNUMBER;
      ADIGATOR.PREOPCOUNT    = NumInVars+1;
      ADIGATOR.VARINFO       = FunctionInfo(FunID).VARINFO;
      ADIGATOR.VARINFO.COUNT = NumInVars+1;
      ADIGATOR.IFDATA        = FunctionInfo(FunID).IFDATA;
      ADIGATORFORDATA        = FunctionInfo(FunID).FORDATA;
      ADIGATOR.FORINFO.COUNT    = 0;  ADIGATOR.FORINFO.EMBEDDEDCOUNT = 0;
      ADIGATOR.FORINFO.INNERLOC = 0;  ADIGATOR.FORINFO.OUTERLOC = 0;
      ADIGATOR.FORINFO.FLAG     = 0;
      ADIGATOR.EMPTYFLAG        = 0;
      ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;
      
      
      ADIGATORDATA            = FunctionInfo(FunID).DATA;
      ADIGATORVARIABLESTORAGE = FunctionInfo(FunID).VARSTORAGE;
      
      if FunID > 1
        FunctionInfo(FunID).Iteration.CallerID = [ADIGATOR.PRINT.FID ADIGATOR.RUNFLAG];
        ADIGATOR.RUNFLAG = 2;
        TempFileName = sprintf('adigatortempsubfunc%1.0f.m',FunID);
        if CallCount == 1
          tempfid = fopen(TempFileName,'w+');
          FunctionInfo(FunID).TempFID = tempfid;
          ADIGATOR.PRINT.FID = tempfid;
        else
          ADIGATOR.PRINT.FID = FunctionInfo(FunID).TempFID;
        end
      end
      PrintFunctionHeader(FunID,FunctionInfo,InVars,InVarStrs);
      if ADIGATOR.DERNUMBER == 1 && any(ADIGATOR.VARINFO.NAMELOCS(1:NumInVars,2) < 0)
        % One of these inputs is used as a direct input to a higher order
        % derivative file - need to mess with names a bit.
        for Icount = 1:NumInVars
          HDfile = ADIGATOR.VARINFO.NAMELOCS(Icount,2);
          if HDfile < 0
            ADIGATOR.VARINFO.NAMELOCS(Icount,2) = 0;
            x = InVars{Icount};
            oldfunstr = cadafuncname(Icount);
            xderiv = x.deriv;
            for Vcount = 1:ADIGATOR.NVAROFDIFF
              if ~isempty(xderiv(Vcount).nzlocs)
                oldderstr = cadadername(funcname,Vcount,Icount);
                fprintf(ADIGATOR.PRINT.FID,[xderiv(Vcount).name,' = ',oldderstr,';\n']);
              end
            end
            fprintf(ADIGATOR.PRINT.FID,[x.func.name,' = ',oldfunstr,';\n']);
          end
        end
      end
      Outputs = Inputs;
      flag = 0;
    else
      % ---------------------- Don't Have to Run ------------------------ %
      FunctionInfo(FunID).Iteration.CallerID = [ADIGATOR.PRINT.FID ADIGATOR.RUNFLAG];
      ADIGATOR.RUNFLAG = 1;
      ADIGATOR.VARINFO = FunctionInfo(FunID).VARINFO;
      ADIGATORFORDATA  = FunctionInfo(FunID).FORDATA;
      ADIGATOR.IFDATA  = FunctionInfo(FunID).IFDATA;
      ADIGATOR.FORINFO = FunctionInfo(FunID).FORINFO;
      ADIGATORVARIABLESTORAGE = FunctionInfo(FunID).VARSTORAGE;
      IterCount = SameFlag;
      FunctionInfo(FunID).Iteration.IterID(CallCount) = IterCount;
      Outputs = FunctionInfo(FunID).Output.Vars(:,IterCount);
      OldLastOcc = FunctionInfo(FunID).VARINFO.LASTOCC;
      [FunctionInfo,Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
      FunctionInfo(FunID).VARINFO.LASTOCC = OldLastOcc;
      flag = 1;
    end
  end
elseif ADIGATOR.RUNFLAG == 1
  %% ~~~~~~~~~~~~~~~~~~~~~~~~~ OVERMAP RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  FunctionInfo = storeOldGlobalParams(FunctionInfo);
  
  if ADIGATOR.EMPTYFLAG
    % ------------------------ Empty Run -------------------------------- %
    Outputs = FunctionInfo(FunID).Output.EmptyVars;
    ADIGATOR.FILE.FUNID    = FunID;
    ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
    [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
    flag = 1;
  elseif ~FunctionInfo(FunID).Iteration.CallCount
    %                          First Call                                 %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    FunctionInfo(FunID).Iteration.CallCount = 1;
    FunctionInfo(FunID).Iteration.IterCount = 1;
    FunctionInfo(FunID).Iteration.IterID    = 1;
    FunctionInfo(FunID).Iteration.CallerID  = GetCallerID(0,FunctionInfo);
    FunctionInfo(FunID).Input.StrucVars     = InVars;
    ADIGATOR.FILE.FUNID    = FunID;
    ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
    ADIGATOR.DERNUMBER     = FunctionInfo(FunID).DERNUMBER;
    ADIGATOR.PREOPCOUNT    = NumInVars+1;
    ADIGATOR.VARINFO       = FunctionInfo(FunID).VARINFO;
    ADIGATOR.IFDATA        = FunctionInfo(FunID).IFDATA;
    ADIGATORFORDATA        = FunctionInfo(FunID).FORDATA;
    ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;
    if ~FunctionInfo(FunID).FunAsLoopFlag
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 0;
      ADIGATOR.FORINFO.INNERLOC = 0;  ADIGATOR.FORINFO.OUTERLOC = 0;
      ADIGATOR.FORINFO.FLAG = 0;
    else
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 1;
      ADIGATOR.FORINFO.OUTERLOC = 1;  ADIGATOR.FORINFO.FLAG = 1;
      ADIGATOR.FORINFO.INNERLOC = 1;
      ADIGATORFORDATA(1).FOR(1).LENGTHS = 1;
      adigatorForIterStart(1,1);
    end
    ADIGATOR.VARINFO.COUNT = NumInVars+1;
    ADIGATOR.EMPTYFLAG             = 0;
    ADIGATORVARIABLESTORAGE = FunctionInfo(FunID).VARSTORAGE;
    for Icount = 1:NumInVars
      cadaOverMap(InVars{Icount});
    end
    Outputs = Inputs;
    flag = 0;
  else
    %                       Second or Higher Call                         %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    if ~isequal(FunctionInfo(FunID).Input.StrucNames,InVarStrs)
      % Check and make sure the input structure is the same
      error(['Currently if using multiple calls to a function, the',...
        ' function must always have the same I/O structure, with the',...
        ' same Structure/Cell patterns of the inputs/outputs.'])
    end
    
    CallCount = FunctionInfo(FunID).Iteration.CallCount + 1;
    FunctionInfo(FunID).Iteration.CallCount = CallCount;
    CallerID  = GetCallerID(0,FunctionInfo);
    FunctionInfo(FunID).Iteration.CallerID(1:length(CallerID),CallCount) = CallerID;
    ADIGATOR.FILE.FUNID    = FunID;
    ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
    % ---------------------- Compare Inputs ----------------------------- %
    for ITcount = 1:FunctionInfo(FunID).Iteration.IterCount
      SameFlag = ITcount;
      for Icount = 1:NumInVars
        InVar1 = InVars{Icount};
        InVar2 = FunctionInfo(FunID).Input.StrucVars{Icount,ITcount};
        InVar1.func.name =[]; InVar2.func.name = [];
        if ~isequal(InVar1.func,InVar2.func)
          SameFlag = 0; break
        else
          for Vcount = 1:ADIGATOR.NVAROFDIFF
            if ~isequal(InVar1.deriv(Vcount).nzlocs,...
                InVar2.deriv(Vcount).nzlocs)
              SameFlag = 0; break
            end
          end
        end
      end
      if SameFlag; break; end
    end
    if SameFlag
      % ------------------- Found Match - Dont Run ---------------------- %
      ADIGATOR.VARINFO = FunctionInfo(FunID).VARINFO;
      ADIGATORFORDATA  = FunctionInfo(FunID).FORDATA;
      ADIGATOR.IFDATA  = FunctionInfo(FunID).IFDATA;
      ADIGATOR.FORINFO = FunctionInfo(FunID).FORINFO;
      ADIGATORVARIABLESTORAGE = FunctionInfo(FunID).VARSTORAGE;
      IterCount = SameFlag;
      FunctionInfo(FunID).Iteration.IterID(CallCount) = IterCount;
      Outputs = FunctionInfo(FunID).Output.Vars(:,IterCount);
      OldLastOcc = FunctionInfo(FunID).VARINFO.LASTOCC;
      [FunctionInfo,Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
      FunctionInfo(FunID).VARINFO.LASTOCC = OldLastOcc;
      flag = 1;
    else
      % ------------------- No Match - Run Iteration -------------------- %
      IterCount = FunctionInfo(FunID).Iteration.IterCount + 1;
      FunctionInfo(FunID).Iteration.IterCount = IterCount;
      FunctionInfo(FunID).Iteration.IterID(CallCount)  = IterCount;
      FunctionInfo(FunID).Input.StrucVars(:,IterCount) = InVars;
      
      ADIGATOR.DERNUMBER     = FunctionInfo(FunID).DERNUMBER;
      ADIGATOR.PREOPCOUNT    = NumInVars+1;
      ADIGATOR.VARINFO       = FunctionInfo(FunID).VARINFO;
      ADIGATOR.IFDATA        = FunctionInfo(FunID).IFDATA;
      ADIGATORFORDATA        = FunctionInfo(FunID).FORDATA;
      ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;
      if FunID == 1
        ADIGATOR.FORINFO.COUNT    = 0;  ADIGATOR.FORINFO.EMBEDDEDCOUNT = 0;
        ADIGATOR.FORINFO.INNERLOC = 0;  ADIGATOR.FORINFO.OUTERLOC = 0;
        ADIGATOR.FORINFO.FLAG = 0;
      else
        ADIGATOR.FORINFO.COUNT    = 1;  ADIGATOR.FORINFO.EMBEDDEDCOUNT = 1;
        ADIGATOR.FORINFO.OUTERLOC = 1;  ADIGATOR.FORINFO.FLAG = 1;
        ADIGATOR.FORINFO.INNERLOC = 1;
        ADIGATORFORDATA(1).FOR(1).LENGTHS = IterCount;
        adigatorForIterStart(1,IterCount);
      end
      ADIGATOR.VARINFO.COUNT = NumInVars+1;
      ADIGATOR.EMPTYFLAG             = 0;
      
      ADIGATORVARIABLESTORAGE = FunctionInfo(FunID).VARSTORAGE;
      for Icount = 1:NumInVars
        cadaOverMap(InVars{Icount});
      end
      Outputs = Inputs;
      flag = 0;
    end
  end
  if ~flag
    ADIGATOR.PREOPCOUNT    = NumInVars+1;
    ADIGATOR.VARINFO.COUNT = NumInVars + 1;
  end
else
  %% ~~~~~~~~~~~~~~~~~~~~~~~~ PRINTING RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  CallCount = FunctionInfo(FunID).Iteration.CallCount;
  if ~CallCount
    %                     Printing this Function                          %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    % ------------------- Set Global Parameters ------------------------- %
    ADIGATOR.VARINFO.LASTOCC = [zeros(size(ADIGATOR.VARINFO.LASTOCC))...
      ADIGATOR.VARINFO.LASTOCC];
    
    ADIGATOR.PRINT.INDENT    = [];
    ADIGATOR.PRINT.FLAG      = 1;
    ADIGATOR.EMPTYFLAG       = 0;
    ADIGATOR.PREOPCOUNT      = NumInVars+1;
    ADIGATOR.FILE.FUNID      = FunID;
    ADIGATOR.FILE.PARENTID   = [ADIGATOR.FILE.PARENTID FunID];
    ADIGATORFORDATA          = FunctionInfo(FunID).FORDATA;
    ADIGATOR.IFDATA          = FunctionInfo(FunID).IFDATA;
    ADIGATORVARIABLESTORAGE  = FunctionInfo(FunID).VARSTORAGE;
    ADIGATOR.STRUCASGN = FunctionInfo(FunID).STRUCASGN;

    
    ADIGATORDATA = FunctionInfo(FunID).DATA;
    if ~FunctionInfo(FunID).FunAsLoopFlag
      ADIGATOR.FORINFO.FLAG          = 0;
      ADIGATOR.FORINFO.INNERLOC      = 0;
      ADIGATOR.FORINFO.OUTERLOC      = 0;
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 0;
    else
      ADIGATOR.FORINFO.FLAG          = 1;
      ADIGATOR.FORINFO.INNERLOC      = 1;
      ADIGATOR.FORINFO.OUTERLOC      = 1;
      ADIGATOR.FORINFO.EMBEDDEDCOUNT = 1;
      adigatorForIterStart(1,1);
    end
    ADIGATOR.VARINFO.COUNT   = NumInVars+1;
    % ------------------- Print Function Header ------------------------- %
    PrintFunctionHeader(FunID,FunctionInfo,InVars,InVarStrs);
    if ADIGATOR.DERNUMBER == 1 && any(ADIGATOR.VARINFO.NAMELOCS(1:NumInVars,2) < 0)
      % One of these inputs is used as a direct input to a higher order
      % derivative file - need to mess with names a bit.
      for Icount = 1:NumInVars
        HDfile = ADIGATOR.VARINFO.NAMELOCS(Icount,2);
        if HDfile < 0
          ADIGATOR.VARINFO.NAMELOCS(Icount,2) = 0;
          x = InVars{Icount};
          oldfunstr = cadafuncname(Icount);
          xderiv = x.deriv;
          for Vcount = 1:ADIGATOR.NVAROFDIFF
            if ~isempty(xderiv(Vcount).nzlocs)
              oldderstr = cadadername(funcname,Vcount,Icount);
              fprintf(ADIGATOR.PRINT.FID,[xderiv(Vcount).name,' = ',oldderstr,';\n']);
            end
          end
          fprintf(ADIGATOR.PRINT.FID,[x.func.name,' = ',oldfunstr,';\n']);
        end
      end
    end
    Outputs = Inputs;
    flag    = 0;
  elseif ADIGATOR.EMPTYFLAG
    % ------------------------ Empty Run -------------------------------- %
    Outputs = FunctionInfo(FunID).Output.EmptyVars;
    ADIGATOR.FILE.FUNID    = FunID;
    ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
    [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
    flag = 1;
  else
    %                     Printing a Parent Function                      %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    fid    = ADIGATOR.PRINT.FID;
    indent = ADIGATOR.PRINT.INDENT;
    % ------------------- Find What Iteration This is --------------------%
    CallerID = [ADIGATOR.FILE.FUNID; ADIGATOR.VARINFO.COUNT];
    for Ccount = 1:size(FunctionInfo(FunID).Iteration.CallerID,2)
      if isequal(CallerID,FunctionInfo(FunID).Iteration.CallerID([1 2],Ccount))
        CallCount = Ccount;
        IterCount = FunctionInfo(FunID).Iteration.IterID(CallCount);
        FunctionInfo(FunID).Iteration.IterCount = IterCount;
        break
      end
    end
    % ------------------------- OverMap Inputs -------------------------- %
    if FunctionInfo(FunID).FunAsLoopFlag
      for Icount = 1:NumInVars
        if InOps(Icount)
          InVar   = InVars{Icount};
          OverLoc = FunctionInfo(FunID).VARINFO.OVERMAP.FOR(Icount,1);
          OverVar = FunctionInfo(FunID).VARSTORAGE.OVERMAP{OverLoc};
          cadaPrintReMap(InVar,OverVar,InOps(Icount));
        end
      end
    end
    
    % ----------------------- Build Call String ------------------------- %
    InCallStr = cell(1,NUMinputs);
    for Icount = 1:NUMinputs
      InCallStr{Icount} = sprintf('cadainput%1.0d,',Icount);
    end
    InCallStr = cell2mat(InCallStr);
    % Check to see if this function is Dependent upon Iterations.
    if FunctionInfo(FunID).Iteration.DepFlag
      % Is dependent - need to add another input so the file knows what
      % iteration this is.
      IterVarStr   = sprintf('adigator%1.0diterID',ADIGATOR.DERNUMBER);
      InCallStr = ['(',InCallStr,IterVarStr,')'];
      % Need to See what all we print out to the IterID
      IterIDStr    = GetCallerID(FunID,FunctionInfo);
      % Print out the IterIDStr.
      fprintf(fid,[indent,IterVarStr,' = ',IterIDStr,';\n']);
    else
      % Is not dependent - do not need to add another input
      InCallStr = ['(',InCallStr(1:end-1),')'];
    end
    NUMoutputs = length(FunctionInfo(FunID).Output.Names);
    if NUMoutputs == 1
      OutCallStr = 'cadaoutput1';
    else
      OutCallStr = cell(1,NUMoutputs);
      for Ocount = 1:NUMoutputs
        OutCallStr{Ocount} = sprintf('cadaoutput%1.0d,',Ocount);
      end
      OutCallStr = cell2mat(OutCallStr);
      OutCallStr = ['[',OutCallStr(1:end-1),']'];
    end
    if FunctionInfo(FunID).DERNUMBER > 1
      FunctionStr = FunctionInfo(FunID).File.Name;
      if length(FunctionStr) > 8 && strcmp(FunctionStr(1:7),'ADiGator')
        FunctionStr = ['ADiGator_',FunctionStr(8:end)];
      end
    else
      FunctionStr = ['ADiGator_',FunctionInfo(FunID).File.Name];
    end
    fprintf(fid,[indent,OutCallStr,' = ',FunctionStr,InCallStr,';\n']);
    if ADIGATOR.OPTIONS.COMMENTS
      fprintf(fid,[indent,'%% Call to function: ',FunctionInfo(FunID).File.Name,'\n']);
    end
    
    % ---------------------- Give Outputs ------------------------------- %
    FunctionInfo = storeOldGlobalParams(FunctionInfo);
    ADIGATOR.FILE.FUNID    = FunID;
    ADIGATOR.FILE.PARENTID = [ADIGATOR.FILE.PARENTID FunID];
    Outputs = FunctionInfo(FunID).Output.Vars(:,IterCount);
    [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs);
    flag = 1;
  end
end

if FunctionInfo(FunID).DERNUMBER > 1
  %% ~~~~~~ LOAD IN PREVIOUS DERIV DATA IF 2ND OR HIGHER CALL ~~~~~~~~~~ %%
  if ~flag
    [varargout{1},FunctionInfo] = getOverloadedData(FunctionInfo,FunID,ADIGATOR.NVAROFDIFF);
  else
    varargout = cell(1);
  end
end

if ADIGATOR.OPTIONS.UNROLL
  ADIGATOR.IFINFO.INNERLOC = 0;
end
end

function [CurVar,Vars,VarStrs,InOps,NumInVars] =...
  FindInputVars(CurVar,CurVarStr,Vars,VarStrs,InOps,NumInVars,PDflag)
% Find Input or Output Variables embedded within cells/structures

global ADIGATOR
if isa(CurVar,'cada')
  % CurVar is overloaded
  NumInVars = NumInVars+1;
  if length(Vars) < NumInVars
    Vars{NumInVars*2,1}    = [];
    VarStrs{NumInVars*2,1} = [];
    InOps(NumInVars*2,1)   = 0;
  end

  
  if ADIGATOR.FILE.FUNID
    varID = CurVar.id;
    %ADIGATOR.VARINFO.LASTOCC(varID,1) = ADIGATOR.VARINFO.COUNT-1;
  else
    varID = NumInVars;
  end
  CurVar.id = NumInVars;
  InOps(NumInVars) = varID;

  if ADIGATOR.RUNFLAG == 2
    funcname = cadafuncname(NumInVars);
    CurVar.func.name = funcname;
    for Vcount = 1:ADIGATOR.NVAROFDIFF
      if ~isempty(CurVar.deriv(Vcount).nzlocs)
        CurVar.deriv(Vcount).name = cadadername(funcname,Vcount,NumInVars);
      end
    end
  end
  Vars{NumInVars}    = CurVar;
  VarStrs{NumInVars} = CurVarStr;
  
elseif isnumeric(CurVar)
  % CurVar is Numeric - turn overloaded
  if length(size(CurVar)) > 2
    error('ADiGator only allows for 2-dimensional arrays')
  end
    NumInVars = NumInVars+1;
    if length(Vars) < NumInVars
      Vars{NumInVars*2,1} = [];
      VarStrs{NumInVars*2,1} = [];
      InOps(NumInVars*2,1) = 0;
    end
    if ADIGATOR.OPTIONS.PREALLOCATE
      Vars{NumInVars}    = CurVar;
      VarStrs{NumInVars} = CurVarStr;
      return
    end
    [xMrow,xNcol] = size(CurVar);
    func.size = [xMrow, xNcol];
    
    varID = NumInVars;
    if ADIGATOR.OPTIONS.AUXDATA == 1 && ~ADIGATOR.FILE.FUNID
      func.zerolocs = find(~CurVar(:));
      if length(func.zerolocs) == xMrow*xNcol
        func.zerolocs = [];
        func.value    = zeros(xMrow,xNcol);
      end
    else
      func.value = CurVar;
    end
    if ADIGATOR.RUNFLAG == 2
      func.name = cadafuncname(NumInVars);
    end
    deriv = struct('name',cell(ADIGATOR.NVAROFDIFF,1),...
      'nzlocs',cell(ADIGATOR.NVAROFDIFF,1));
    CurVar = cada(varID,func,deriv);
    
    Vars{NumInVars}    = CurVar;
    VarStrs{NumInVars} = CurVarStr;
elseif isstruct(CurVar)
  % CurVar is Structure - search all elements/fields
  if length(size(CurVar)) > 2
    error('Can only use two dimensional Structure Arrays')
  end
  sNum        = numel(CurVar);
  % Check for pp or adigatorpp2 structures
%   if (isfield(CurVar,'form') && isfield(CurVar,'breaks') && isfield(CurVar,'coefs') && ...
%       isfield(CurVar,'pieces') && isfield(CurVar,'order') && isfield(CurVar,'dim')) || ...
%       (isfield(CurVar,'form') && isfield(CurVar,'xbreaks') && isfield(CurVar,'ybreaks') ...
%       && isfield(CurVar,'coefs') && isfield(CurVar,'xorder') && isfield(CurVar,'yorder'))
%     for Icount = 1:size(CurVar,1)
%       for Jcount = 1:size(CurVar,2)
%         func.size = [1 1];
%         NumInVars = NumInVars+1;
%         if length(Vars) < NumInVars
%           Vars{NumInVars*2,1} = [];
%           VarStrs{NumInVars*2,1} = [];
%           InOps(NumInVars*2,1) = 0;
%         end
%         if ADIGATOR.OPTIONS.PREALLOCATE
%           Vars{NumInVars}    = CurVar;
%           VarStrs{NumInVars} = CurVarStr;
%           return
%         end
%         varID = NumInVars;
%         func.value = CurVar(Icount,Jcount);
%         if ADIGATOR.RUNFLAG == 2
%           func.name = cadafuncname(NumInVars);
%         end
%         func.pp = 1;
%         deriv = struct('name',cell(ADIGATOR.NVAROFDIFF,1),...
%           'nzlocs',cell(ADIGATOR.NVAROFDIFF,1));
%         CurVar = cada(varID,func,deriv);
%         Vars{NumInVars}    = CurVar;
%         if sNum == 1
%           VarStrs{NumInVars} = CurVarStr;
%         else
%           VarStrs{NumInVars} = sprintf([CurVarStr,'(%1.0f,%1.0f)'],Icount,Jcount);
%         end
%       end
%     end
%     return
%   end
  CurVar      = orderfields(CurVar);
  StructNames = fieldnames(CurVar);
  if sNum > 1
    for Icount = 1:size(CurVar,1)
      for Jcount = 1:size(CurVar,2)
        for Scount = 1:length(StructNames)
          NewVar    = CurVar(Icount,Jcount).(StructNames{Scount});
          NewVarStr = sprintf([CurVarStr,'(%1.0f,%1.0f).',StructNames{Scount}],Icount,Jcount);
          [CurVar(Icount,Jcount).(StructNames{Scount}),Vars,VarStrs,InOps,NumInVars] =...
            FindInputVars(NewVar,NewVarStr,Vars,VarStrs,InOps,NumInVars,PDflag);
        end
      end
    end
  else
    for Scount = 1:length(StructNames)
      NewVar    = CurVar.(StructNames{Scount});
      NewVarStr = [CurVarStr,'.',StructNames{Scount}];
      [CurVar.(StructNames{Scount}),Vars,VarStrs,InOps,NumInVars] =...
        FindInputVars(NewVar,NewVarStr,Vars,VarStrs,InOps,NumInVars,PDflag);
    end
  end
elseif iscell(CurVar)
  % CurVar is Cell - search all elements
  cSize = size(CurVar);
  if length(cSize) > 2
    error('ADiGator only allows for 2-dimensional cell arrays');
  end
  cRow = cSize(1); cCol = cSize(2);
    for Icount = 1:cRow
      for Jcount = 1:cCol
        NewVar = CurVar{Icount,Jcount};
        NewVarStr = sprintf([CurVarStr,'{%1.0d,%1.0d}'],Icount,Jcount);
        [CurVar{Icount,Jcount},Vars,VarStrs,InOps,NumInVars] =...
          FindInputVars(NewVar,NewVarStr,Vars,VarStrs,InOps,NumInVars,PDflag);
      end
    end
else
  %error('Variable inputs/outputs may only be numeric or cell/structure arrays of numeric values.')
end
end

function FunctionInfo = storeOldGlobalParams(FunctionInfo)
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE ADIGATORDATA
if ~isempty(ADIGATOR.FILE.PARENTID)
  ParentID = ADIGATOR.FILE.PARENTID(end);
  FunctionInfo(ParentID).VARINFO    = ADIGATOR.VARINFO;
  FunctionInfo(ParentID).FORINFO    = ADIGATOR.FORINFO;
  FunctionInfo(ParentID).FORDATA    = ADIGATORFORDATA;
  FunctionInfo(ParentID).IFDATA     = ADIGATOR.IFDATA;
  FunctionInfo(ParentID).VARSTORAGE = ADIGATORVARIABLESTORAGE;
  if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.RUNFLAG
    FunctionInfo(ParentID).DATA = ADIGATORDATA;
  end
end
end

function CallerID = GetCallerID(CFunID,FunctionInfo)
% Find the Unique CALLERID from the function.
global ADIGATOR ADIGATORFORDATA

if ~CFunID
  % This gets the CallerID for the Pre-Print Run
  
  % First two things that define the Function Call are the Caller
  % Functions FunID  and the Operation Count in the Caller Function.
  CallerID = [ADIGATOR.FILE.FUNID; ADIGATOR.VARINFO.COUNT];

  if ADIGATOR.FORINFO.FLAG
    % If we are in a FOR loop, we also need the Iteration Count of every
    % FOR loop this call is nested in to Uniquely define the call. - If
    % called from another sub-function, the outermost loop will be that
    % function.
    OFloc = ADIGATOR.FORINFO.OUTERLOC; 
    IFloc = ADIGATOR.FORINFO.INNERLOC;
    for Fcount = 1:length(ADIGATORFORDATA(OFloc).FOR)
      FORLOCS = ADIGATORFORDATA(OFloc).FOR(Fcount).LOCS;
      if FORLOCS(1,end) == IFloc
        % Found the set of FORLOCS that describes the nested loop we are in
        break
      end
    end
    CallerID(end+size(FORLOCS,2),1) = 0;
    
    for Fcount = 1:size(FORLOCS,2)
      % Get the Iteration Count of each Nested FOR loop.
      CallerID(Fcount+2) = ADIGATORFORDATA(FORLOCS(1,Fcount)).COUNT.ITERATION;
    end
  end
else
  % This looks at the CallerID whenever a function is dependent and needs
  % to print out an ID check in the Printing Evaluation case where the
  % sub function is being called by another function.
  
  CALLERID  = FunctionInfo(CFunID).Iteration.CallerID;
  NumID     = size(CALLERID,1);
  CallerDep = zeros(size(CALLERID,1),1);
  for Icount = 1:NumID
    if sum(CALLERID(Icount,1) ~= CALLERID(Icount,:))
      CallerDep(Icount) = 1;
    end
  end
  CallerDep = logical(CallerDep);
  
  % Get the String which contains the dependent rows of CALLERID in this
  % calling function.
  CallerFunIDs  = unique(CALLERID(1,:));
  ThisCallerLoc = CallerFunIDs == ADIGATOR.FILE.FUNID;
  if ~isfield(FunctionInfo(CFunID),'CallerIndices') || isempty(FunctionInfo(CFunID).CallerIndices)
    FunctionInfo(CFunID).CallerIndices = cell(length(CallerFunIDs),2);
  end

  CallerIDstr = FunctionInfo(CFunID).CallerIndices{ThisCallerLoc,1};
  IterIDstr   = FunctionInfo(CFunID).CallerIndices{ThisCallerLoc,2};
  if isempty(CallerIDstr)
    CallerIDstr = cadaindprint(CALLERID(CallerDep,:));
    FunctionInfo(CFunID).CallerIndices{ThisCallerLoc,1} = CallerIDstr;
    IterIDstr = cadaindprint(FunctionInfo(CFunID).Iteration.IterID);
    FunctionInfo(CFunID).CallerIndices{ThisCallerLoc,2} = IterIDstr;
  end
  
  % Build the string that we will use to determine the iteration count.
  NumDep  = sum(CallerDep);
  CallerID = cell(1,NumDep); DepCount = 0;
  % Caller Function ID
  if CallerDep(1)
    DepCount = DepCount+1;
    CallerID{1} =  sprintf(['%1.0d == ',CallerIDstr,'(1,:) & ']);
  end
  % Caller Function Operation Count
  if CallerDep(2)
    DepCount = DepCount+1;
    CallerID{DepCount} = sprintf(['%1.0d == ',CallerIDstr,...
      '(%1.0d,:) & '],ADIGATOR.VARINFO.COUNT,DepCount);
  end
  % Caller Function FOR counts
  if ADIGATOR.FORINFO.FLAG
    OFloc = ADIGATOR.FORINFO.OUTERLOC; 
    IFloc = ADIGATOR.FORINFO.INNERLOC;
    for Fcount = 1:length(ADIGATORFORDATA(OFloc).FOR)
      FORLOCS = ADIGATORFORDATA(OFloc).FOR(Fcount).LOCS;
      if FORLOCS(1,end) == IFloc
        % Found the set of FORLOCS that describes the nested loop we are in
        break
      end
    end
  end
  for Ccount = 3:NumID
    if CallerDep(Ccount)
      DepCount = DepCount+1;
      if ADIGATOR.FORINFO.FLAG
        % This call is from within a FOR loop.
        if Ccount == 3 && ADIGATOR.FILE.FUNID > 1 && ...
            ~any(FunctionInfo(ADIGATOR.FILE.FUNID).Iteration.DepFlag)
          % The calling function only has one iteration.
          CallerID{DepCount} = sprintf(['1 == ',CallerIDstr,...
            '(%1.0d,:) & '],DepCount);
        else
          CallerID{DepCount} = sprintf([...
            ADIGATORFORDATA(FORLOCS(1,Ccount-2)).COUNTNAME,' == ',...
            CallerIDstr,'(%1.0d,:) & '],DepCount);
        end
      else
        % This particular call isnt from within a FOR loop.
        CallerID{DepCount} = sprintf(['0 == ',CallerIDstr,...
            '(%1.0d,:) & '],DepCount);
      end
    end
  end
  CallerID = cell2mat(CallerID);
  CallerID = [IterIDstr,'(',CallerID(1:end-3),')'];
end
return
end

function PrintFunctionHeader(FunID,FunctionInfo,InVars,InVarStrs)
global ADIGATOR ADIGATORFORDATA

fid = ADIGATOR.PRINT.FID;
NumInVars = length(InVars);
InNames  = FunctionInfo(FunID).Input.Names;
OutNames = FunctionInfo(FunID).Output.Names;
NUMout   = length(OutNames);
NUMinputs = length(InNames);

if NUMout == 1
  OutStr = OutNames{1};
else
  OutStr = cell(1,NUMout);
  for Ocount = 1:NUMout
    OutStr{Ocount} = [OutNames{Ocount},','];
  end
  OutStr = cell2mat(OutStr);
  OutStr = ['[',OutStr(1:end-1),']'];
end
InStr = cell(1,NUMinputs);
for Icount = 1:NUMinputs
  InStr{Icount} = [InNames{Icount},','];
end
InStr = cell2mat(InStr);
MatFileName = ['ADiGator_',ADIGATOR.PRINT.FILENAME];
if FunID == 1
  % Main Function
  FileName    = ADIGATOR.PRINT.FILENAME;
  fprintf(fid,['function ',OutStr,' = ',FileName,'(',InStr(1:end-1),')\n']);
  fprintf(fid,['global ',MatFileName,'\n']);
  fprintf(fid,['if isempty(',MatFileName,'); ADiGator_LoadData(); end\n']);
else
  if FunctionInfo(FunID).DERNUMBER > 1
    FileName = FunctionInfo(FunID).File.Name;
    if length(FileName) > 8 && strcmp(FileName(1:7),'ADiGator')
      FileName = ['ADiGator_',FileName(8:end)];
    end
  elseif ADIGATOR.OPTIONS.UNROLL && FunID > 1
    IterCount = FunctionInfo(FunID).Iteration.IterCount;
    FileName = sprintf(['ADiGator_',FunctionInfo(FunID).File.Name,'%1.0f'],IterCount);
  else
    FileName = ['ADiGator_',FunctionInfo(FunID).File.Name];
  end
  if ~ADIGATOR.OPTIONS.UNROLL && FunctionInfo(FunID).Iteration.DepFlag
    fprintf(fid,['function ',OutStr,' = ',FileName,'(',InStr,...
      ADIGATORFORDATA(1).COUNTNAME,')\n']);
  else
    fprintf(fid,['function ',OutStr,' = ',FileName,'(',InStr(1:end-1),')\n']);
  end
  fprintf(fid,['global ',MatFileName,'\n']);
end
% Print Data References
for Dcount = 1:ADIGATOR.DERNUMBER
  fprintf(fid,['Gator%1.0dIndices = ',MatFileName,'.',FileName,...
    '.Gator%1.0dIndices;\n'],Dcount,Dcount);
  fprintf(fid,['Gator%1.0dData = ',MatFileName,'.',FileName,...
    '.Gator%1.0dData;\n'],Dcount,Dcount);
end


fprintf(fid,'%% ADiGator Start Derivative Computations\n');
if ~ADIGATOR.OPTIONS.UNROLL && length(FunctionInfo(FunID).Iteration.DepFlag) == 2 &&...
    ~FunctionInfo(FunID).Iteration.DepFlag(2)
  fprintf(fid,'adigator%1.0dfuncount = adigator%1.0dfuncount;\n',...
    ADIGATOR.DERNUMBER,ADIGATOR.DERNUMBER-1);
end


for Icount = 1:NumInVars
  cadaOverMap(InVars{Icount});
end

return
end

function [OverLoadedData,FunctionInfo] = getOverloadedData(FunctionInfo,FunID,NUMvod)
PreviousDerivData = FunctionInfo(FunID).PreviousDerivData;
NumPrevDerivs     = length(PreviousDerivData.Derivative);
if ~isfield(PreviousDerivData,'OverLoadedData')
  TempDeriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  for Dcount = 1:NumPrevDerivs
    DIndexStr = sprintf('Gator%1.0dIndices',Dcount);
    iIndexFields = fieldnames(PreviousDerivData.(DIndexStr));
    for Icount = 1:length(iIndexFields)
      IndexStr = iIndexFields{Icount};
      TempFunc.name = [DIndexStr,'.',IndexStr];
      TempFunc.value = PreviousDerivData.(DIndexStr).(IndexStr);
      TempFunc.size = size(TempFunc.value);
      OverLoadedData.(DIndexStr).(IndexStr) = cada([],TempFunc,TempDeriv);
    end
    
    DDataStr = sprintf('Gator%1.0dData',Dcount);
    IDataFields = fieldnames(PreviousDerivData.(DDataStr));
    for Icount = 1:length(IDataFields)
      DataStr = IDataFields{Icount};
      TempFunc.name = [DDataStr,'.',DataStr];
      TempFunc.value = PreviousDerivData.(DDataStr).(DataStr);
      if (isstruct(TempFunc.value) && isfield(TempFunc.value,'form') && ...
          ischar(TempFunc.value.form) && strcmp(TempFunc.value.form,'adigatorpp2'))||...
        (isstruct(TempFunc.value) && isfield(TempFunc.value,'pp') && ...
          ischar(TempFunc.value.form) && strcmp(TempFunc.value.form,'pp'))
        TempFunc.pp = 1;
      end
      TempFunc.size = size(TempFunc.value);
      OverLoadedData.(DDataStr).(DataStr) = cada([],TempFunc,TempDeriv);
    end
  end
  FunctionInfo(FunID).PreviousDerivData.OverLoadedData = OverLoadedData;
else
  OverLoadedData = PreviousDerivData.OverLoadedData;
end
end