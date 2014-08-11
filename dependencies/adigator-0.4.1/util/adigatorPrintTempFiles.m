function [ForCount, IfCount] = adigatorPrintTempFiles(Ffid,Tfid,FlowInfo,...
  DerNumber,ForCount,FunStrChecks)
% adigatorPrintTempFiles(Ffid,FlowInfo)
% This function is used to make blocks of temporary functions, which
% contain lines of user code interspersed with statements which call
% adigatorVarAnalyzer in order to read and print out the derivative function 
% properly.
% This routine is called from adigator.m and calls no routines which are not
% sub-routines of adigatorPrintTempFiles.m itself.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

% ----------------------------------------------------------------------- %
%                                MAIN                                     %
% ----------------------------------------------------------------------- %

IfCount  = 0;
indent   = [];
MainFlag = ~ForCount;
ForIfFlag  = [0 0 0];
StartLocation = FlowInfo.StartLocation;
if ~isempty(FlowInfo.Children)
  % Print the first block of code.
  EndLocation = FlowInfo.Children(1).StartLocation;
  PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
  % Print the children blocks of code.
  [StartLocation,ForCount,IfCount,~] = FlowPrint(Ffid,FlowInfo.Children,...
    Tfid,indent,ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag);
end
EndLocation = FlowInfo.EndLocation;
% Print the last block of code.
PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
return
end

function [EndLocation,ForCount,IfCount,ForIfFlag] = FlowPrint(Ffid,FlowInfo,Tfid,...
  indent,ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag)
global ADIGATOR
% Function Recursively calls itself to get through the FlowInfo structure
% and print out all of the temporary files.
FlowSize = length(FlowInfo);
LastBro  = 0;
for FlowCount = 1:FlowSize
  if strcmp(FlowInfo(FlowCount).Type,'if')
    % ------------------------------------------------------------------- %
    %                         IF STATEMENTS - START                       %
    % ------------------------------------------------------------------- %
    IfCount       = IfCount+1;
    ForIfFlag(2)  = ForIfFlag(2)+1;
    BroCount = 1;
    CurIfCount = IfCount;
    if FlowCount ~= 1
      % Need to get the section between previous Flow Control Statment and
      % this one.
      StartLocation = EndLocation; % Starting at the End of the last section
      EndLocation = FlowInfo(FlowCount).StartLocation;
      PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    end
    fprintf(Tfid,[indent,'%% ADiGator IF Statement #%1.0d: START\n'],IfCount);
    
    % --------------- Print Out the Conditional Variables --------------- %
    CondVarStr = getIfForStatement(Ffid,FlowInfo(FlowCount).StartLocation);
    CadaVarStr = 'cadaconditional1';
    if DerNumber == 1
      fprintf(Tfid,[indent,CadaVarStr,' = ',CondVarStr,';\n']);
      CondVarStr = FindDoinkers(CondVarStr);
      fprintf(Tfid,[indent,CadaVarStr,' = adigatorVarAnalyzer(''',CadaVarStr,...
        ' = ',CondVarStr,''',',CadaVarStr,',''',CadaVarStr,''',0);\n']);
    end
    CondVarStrs = cell(1,FlowSize-FlowCount);
    CondVarStrs{1} = [CadaVarStr,','];
    % Get Any ElseIf variables
    for FlowCount2 = FlowCount+1:FlowSize
      if strcmp(FlowInfo(FlowCount2).Type,'elseif')
        BroCount = BroCount+1;
        if DerNumber == 1
          CondVarStr = getIfForStatement(Ffid,FlowInfo(FlowCount2).StartLocation);
          CadaVarStr = sprintf('cadaconditional%1.0d',BroCount);
          fprintf(Tfid,[indent,CadaVarStr,' = ',CondVarStr,';\n']);
          CondVarStr = FindDoinkers(CondVarStr);
          fprintf(Tfid,[indent,CadaVarStr,' = adigatorVarAnalyzer(''',CadaVarStr,...
            ' = ',CondVarStr,''',',CadaVarStr,',''',CadaVarStr,''',0);\n']);
        end
        CondVarStrs{BroCount} = [CadaVarStr,','];
      elseif strcmp(FlowInfo(FlowCount2).Type,'else')
        BroCount = BroCount+1;
        CondVarStrs{BroCount} = '[],';
      else
        break
      end
    end
    LastBro = BroCount;
    % ---------------- Print Out the IF Initialize Statement ------------ %
    CondVarStrs = CondVarStrs(1:BroCount);
    CondVarStrs = cell2mat(CondVarStrs);
    if ADIGATOR.OPTIONS.UNROLL && (ForIfFlag(2) == 1 || ForIfFlag(3))
        fprintf(Tfid,[indent,'for adigatorIfPrint%1.0f = adigatorIfLooper(%1.0f)\n'],...
          IfCount,IfCount);
        fprintf(Tfid,[indent,'adigatorIfLooperi(adigatorIfPrint%1.0f,%1.0f);\n'],IfCount,IfCount);
    end
    fprintf(Tfid,[indent,'adigatorIfInitialize(%1.0f,',...
      CondVarStrs(1:end-1),');\n'],IfCount);
    if ADIGATOR.OPTIONS.UNROLL
      fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
        'adigatorIfIterStart(%1.0f,1);%%#ok<NASGU>\n'],IfCount);
      fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr);',...
        ' cellfun(@eval,adigatorIfEvalStr); end\n']);
    else
      fprintf(Tfid,[indent,'adigatorIfIterStart(%1.0f,1);\n'],IfCount);
    end
    BroCount = 1;
  elseif strcmp(FlowInfo(FlowCount).Type,'elseif') || strcmp(FlowInfo(FlowCount).Type,'else')
    % ------------------------------------------------------------------- %
    %                  ELSEIF/ELSE STATEMENTS  - START                    %
    % ------------------------------------------------------------------- %
    BroCount = BroCount+1;
    % ------------------ Print Out the ELSEIF Statement ----------------- %
        fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
      'adigatorIfIterStart(%1.0f,%1.0f);%%#ok<NASGU>\n'],CurIfCount,BroCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr);',...
      ' cellfun(@eval,adigatorIfEvalStr); end\n']);
  elseif strcmp(FlowInfo(FlowCount).Type,'for')
    % ------------------------------------------------------------------- %
    %                      FOR STATEMENTS  - START                        %
    % ------------------------------------------------------------------- %
    ForCount      = ForCount+1;
    ForIfFlag(1)  = ForIfFlag(1)+1;
    if ForIfFlag(2); ForIfFlag(3) = ForIfFlag(3)+1; end
    CurForCount   = ForCount;
    if FlowCount ~= 1
      % Need to get the section between previous Flow Control Statment and
      % this one.
      StartLocation = EndLocation; % Starting at the End of the last section
      EndLocation = FlowInfo(FlowCount).StartLocation;
      PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    end
    fprintf(Tfid,[indent,'%% ADiGator FOR Statement #%1.0d: START\n'],ForCount);
    % -------------------- Get the Loop Variable ------------------------ %
    LoopStr = getIfForStatement(Ffid,FlowInfo(FlowCount).StartLocation);
    % See if we need to print it out or not
    EqLoc = strfind(LoopStr,'=');
    if ~isempty(EqLoc)
      LoopStrLHS = strtrim(LoopStr(1:EqLoc(1)-1));
      LoopStrRHS = strtrim(LoopStr(EqLoc(1)+1:end));
      LoopVar    = sprintf('cadaforvar%1.0d',ForCount);
      LoopCount  = sprintf('cadaforcount%1.0d',ForCount);
      if DerNumber == 1
        LoopVarStr = [LoopVar,' = ',LoopStrRHS,';'];
        fprintf(Tfid,[indent,LoopVarStr,'\n']);
        fprintf(Tfid,[indent,LoopVar,' = adigatorVarAnalyzer(''',FindDoinkers(LoopVarStr),''',',...
          LoopVar,',''',LoopVar,''',0);\n']);
      else
        LoopVar = LoopStrRHS;
        LoopCount = LoopStrLHS;
      end
    else
      errlink = GenErrorLink(Ffid,FlowInfo(FlowCount).StartLocation(1));
      error(['???Unable to parse ',LoopStr,': No Equal Sign at: ',errlink])
    end
    
    % --------------- Print the FOR Initialize Statement ---------------- %
    if ForIfFlag(1) == 1 && ~ADIGATOR.OPTIONS.UNROLL
      % Main function on an outer loop
      fprintf(Tfid,[indent,'[adigatorForVariable%1.0d, adigatorForEvalStr, adigatorForEvalVar]',...
        ' = adigatorForInitialize(%1.0d,',LoopVar,');%%#ok<NASGU>\n'],ForCount,ForCount);
      fprintf(Tfid,[indent,'if ~isempty(adigatorForEvalStr)\n']);
      fprintf(Tfid,[indent,'    cellfun(@eval,adigatorForEvalStr);\n']);
      fprintf(Tfid,[indent,'end\n']);
    else
      fprintf(Tfid,[indent,'adigatorForVariable%1.0d = adigatorForInitialize(%1.0d,',...
        LoopVar,');\n'],ForCount,ForCount);
    end
    fprintf(Tfid,[indent,'for adigatorForVariable%1.0di = adigatorForVariable%1.0d\n'],...
      ForCount,ForCount);
    fprintf(Tfid,[indent,LoopCount,' = adigatorForIterStart(%1.0d,',...
      'adigatorForVariable%1.0di);\n'],ForCount,ForCount);
    if DerNumber == 1
      fprintf(Tfid,[indent,LoopStrLHS,' = ',LoopVar,'(:,',LoopCount,');\n']);
      fprintf(Tfid,[indent,LoopStrLHS,' = adigatorVarAnalyzer(''',LoopStrLHS,...
        ' = ',LoopVar,'(:,',LoopCount,');'',',LoopStrLHS,',''',LoopStrLHS,''',0);\n']);
    end
  end
  
  % --------------------------------------------------------------------- %
  %                      CALCULATIONS WITHIN THE FLOW CONTROL             %
  % --------------------------------------------------------------------- %
  indent = [indent,'    ']; %#ok<AGROW>
  StartLocation = FlowInfo(FlowCount).StartLocation;
  if ~isempty(FlowInfo(FlowCount).Children)
    % Print out the block leading to first child.
    EndLocation = FlowInfo(FlowCount).Children(1).StartLocation;
    PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    % Print the children blocks of code.
    [StartLocation,ForCount,IfCount,ForIfFlag] = ...
      FlowPrint(Ffid,FlowInfo(FlowCount).Children,Tfid,indent,...
      ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag);
  end
  EndLocation = FlowInfo(FlowCount).EndLocation;
  % Print the last block of code for this section.
  PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber)

  % --------------------------------------------------------------------- %
  %                      END OF THE FLOW CONTROL BLOCK                    %
  % --------------------------------------------------------------------- %
  indent(1:4) = [];  
  if strcmp(FlowInfo(FlowCount).Type,'for')
    fprintf(Tfid,[indent,'[adigatorForEvalStr, adigatorForEvalVar]',...
      '= adigatorForIterEnd(%1.0d,adigatorForVariable%1.0di);\n'],CurForCount,CurForCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorForEvalStr)\n']);
    fprintf(Tfid,[indent,'    cellfun(@eval,adigatorForEvalStr);\n']);
    fprintf(Tfid,[indent,'end\n']);
    fprintf(Tfid,[indent,'end\n']);
    fprintf(Tfid,[indent,'%% ADiGator FOR Statement #%1.0d: END\n'],CurForCount);
    ForIfFlag(1) = ForIfFlag(1)-1;
    if ForIfFlag(2); ForIfFlag(3) = ForIfFlag(3) - 1; end
  elseif BroCount == LastBro
    fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
      'adigatorIfIterEnd(%1.0f,%1.0f);%%#ok<NASGU>\n'],CurIfCount,BroCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr);',...
      ' cellfun(@eval,adigatorIfEvalStr); end\n']);
    if ADIGATOR.OPTIONS.UNROLL && (ForIfFlag(2) == 1 || ForIfFlag(3))
      fprintf(Tfid,[indent,'end\n']);
    end
    fprintf(Tfid,[indent,'%% ADiGator IF Statement #%1.0d: END\n'],CurIfCount);
    ForIfFlag(2) = ForIfFlag(2)-1;
  else
    fprintf(Tfid,[indent,'adigatorIfIterEnd(%1.0f,%1.0f);\n'],CurIfCount,BroCount);
  end
end

return
end

function PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber)
% Read from file Ffid, print to temporary function.
% Temporary function is named by the EndLocation(3).
global ADIGATOR
VAstr = 'adigatorVarAnalyzer';
fseek(Ffid,StartLocation(4),-1);

% ----------------First Line If There are Multiple Evaluations------------%
FunStrFULL = fgets(Ffid);
MajorLineCount = StartLocation(1);
MinorLineStart = StartLocation(2)+1;

while MajorLineCount <= EndLocation(1) && ~isnumeric(FunStrFULL)
  FunStrFull = strtrim(FunStrFULL);
  if ~isempty(FunStrFull)
    if length(FunStrFull) > 3 && ~strcmp(FunStrFull(1),'%')
      multlines1  = strcmp(FunStrFull(end-2:end),'...');
      osquarelocs = strfind(FunStrFull,'[');
      csquarelocs = strfind(FunStrFull,']');
      multlines2  = length(osquarelocs) > length(csquarelocs);
      while multlines1 || multlines2
        % Single command spanning multiple lines - look for a comment at
        % end of this line
        commentloc = strfind(FunStrFull,'%');
        if ~isempty(commentloc)
          % Just trash the comment - no good way of keeping it
          FunStrFull = strtrim(FunStrFull(1:commentloc(1)-1));
        end
        if multlines1 && length(FunStrFull) > 6 && strcmp(FunStrFull(1:7),'global ')
          % ... at end of line with global
           FunStrFull = [FunStrFull(1:end-3),' '];
        elseif multlines1
          % ... at end of line
          FunStrFull = FunStrFull(1:end-3);
        elseif strcmp(FunStrFull(end),',')
          % Building a Matrix and have ',' at end for some reason
          FunStrFull = [FunStrFull(1:end-1),';'];
        elseif ~strcmp(FunStrFull(end),';')
          % Building a Matrix and dont have the vertcat operator
          FunStrFull = [FunStrFull,';']; %#ok<AGROW>
        end
        FunStrFull = strtrim([FunStrFull,fgets(Ffid)]);
        MajorLineCount = MajorLineCount + 1;
        multlines1  = strcmp(FunStrFull(end-2:end),'...');
        osquarelocs = strfind(FunStrFull,'[');
        csquarelocs = strfind(FunStrFull,']');
        multlines2  = length(osquarelocs) > length(csquarelocs);
      end
    end
    [FunStr,NUMFunStr] = adigatorSeperateFunLines(FunStrFull);
    if MajorLineCount == EndLocation(1)
      NUMFunStr = EndLocation(2)-1;
    end
    % -----------------Work on Function Lines--------------------------
    for adigatorFScount = MinorLineStart:NUMFunStr
      FunStri = FunStr{adigatorFScount};
      StrLength = length(FunStri);
      FunStri = strtrim(FindComments(FunStri));
      SlashLocs = strfind(FunStri,'\');
      if ~isempty(SlashLocs)
        FunStri = strtrim(FindSlashes(FunStri,SlashLocs));
      end
      EqualLoc = regexp(FunStri,'[^=><~]=[^=]')+1;
      if strcmp(FunStri(1),'%')
        % COMMENT
        FunStri = FindDoinkers(FunStri);
        fprintf(Tfid,[indent,VAstr,'(''',FunStri,''');\n']);
      elseif ~isempty(EqualLoc)
        % Some sort of Assignment/Calculation
        FunStri = CheckFunctionCall(FunStri,FunStrChecks,Tfid,indent,DerNumber,MajorLineCount);
        if isempty(FunStri)
          continue
        end
        % --- Check For Structure Subsref on RHS ---
        svrcount = 0; svrtype = 1;
        FunStrRHS = strtrim(FunStri(EqualLoc+1:end));
        [strucloc1,strucloc2] = regexpi(FunStrRHS,'\([^\)\(]*[a-z][^\)\(]*\)\.[a-z]',...
          'start','end','once');
        strucloc2 = strucloc2-1;
        [cellyloc1,cellyloc2] = regexpi(FunStrRHS,'\{[^\}\{]*[a-z][^\}\{]*\}',...
          'start','end','once');
        if isempty(strucloc1) || (~isempty(cellyloc1) && strucloc1 > cellyloc1)
          strucloc1 = cellyloc1; strucloc2 = cellyloc2; svrtype = 0;
        end
        while ~isempty(strucloc1)
          fprintf(Tfid,[indent,'%% ',FunStrRHS,'\n']);
          svrcount = svrcount+1;
          if svrtype
            %svendloc = regexp(FunStrRHS,'\([^\)]*[a-z].*\)\.\w+','end','once');
            svendloc = regexp(FunStrRHS(strucloc2+1:end),'\W','once','start')+strucloc2-1;
            svstartloc = regexp(FunStrRHS(1:strucloc1-2),'[\w\(\)\{\}\.]+$','start','once');
            StrucRef.i    = FunStrRHS(strucloc1+1:strucloc2-2);
          else
            svendloc = strucloc2;
            svstartloc = regexp(FunStrRHS(1:strucloc1-1),'[\w\(\)\{\}\.]+$','start','once');
            StrucRef.i    = FunStrRHS(strucloc1+1:strucloc2-1);
          end
          ADIGATOR.SVACOUNT = ADIGATOR.SVACOUNT+1;
          if DerNumber == 1
            StrucRef.new  = sprintf('adigatorStrucRef%1.0f',ADIGATOR.SVACOUNT);
          else
            StrucRef.new = strtrim(FunStri(1:EqualLoc-1));
          end
          StrucRef.old  = FunStrRHS(svstartloc:svendloc);
          openloc = strfind(StrucRef.old,'(');
          closeloc = strfind(StrucRef.old,')');
          if length(openloc) > length(closeloc)
            svstartloc = openloc(length(openloc)-length(closeloc))+svstartloc;
            StrucRef.old  = FunStrRHS(svstartloc:svendloc);
          end
          StrucRef.s    = FunStrRHS(svstartloc:strucloc1-1);
          %FunStrRHS  = [StrucRef.new,FunStrRHS(svendloc+1:end)];
          strucloc1 = strfind(FunStrRHS,StrucRef.old);
          FunStrRHS = [FunStrRHS(1:strucloc1-1) StrucRef.new,...
            FunStrRHS(strucloc1+length(StrucRef.old):end)];
          StrucRef.old = FindDoinkers(StrucRef.old);
          fprintf(Tfid,[indent,'[',StrucRef.new,',',StrucRef.s,'] = adigatorCellStrucRef(%1.0f,',...
            StrucRef.s,',''',StrucRef.new,' = ',StrucRef.old,';'',',StrucRef.i,');\n'],...
            ADIGATOR.SVACOUNT);
          svrtype = 1; % 1 for structure, 0 for cell
          [strucloc1,strucloc2] = regexpi(FunStrRHS,'\([^\)\(]*[a-z][^\)\(]*\)\.[a-z]',...
            'start','end','once');
          strucloc2 = strucloc2-1;
          [cellyloc1,cellyloc2] = regexpi(FunStrRHS,'\{[^\}]*[a-z][^\}]*\}',...
            'start','end','once');
          if isempty(strucloc1) || (~isempty(cellyloc1) && strucloc1 > cellyloc1)
            strucloc1 = cellyloc1; strucloc2 = cellyloc2; svrtype = 0;
          end
        end
        if svrcount
          ADIGATOR.OPTIONS.PREALLOCATE = 1;
          if DerNumber == 1
            FunStri = [FunStri(1:EqualLoc),FunStrRHS];
          else
            continue
          end
        end
          
          
        % ---Work on LHS of equal Sign---
        VarStr = strtrim(FunStri(1:EqualLoc-1));
        if strcmp(VarStr(1),'[')
          % Multiple Assignments
          VarStrings = SeperateOutputStrings(strtrim(VarStr(2:end-1)),0);
          NumOutVars = length(VarStrings);
        else
          % Single Assignment
          VarStrings = cell(1); VarStrings{1} = VarStr;
          NumOutVars = 1;
        end
        SubsFlags = zeros(NumOutVars,1);
        
        svacount = ADIGATOR.SVACOUNT;
        for Vcount = 1:NumOutVars
          VarStr = VarStrings{Vcount};
          % Check for Structure Assignments
          svatype = 1; % 1 for structure, 0 for cell
          [strucloc1,strucloc2] = regexpi(VarStr,'\([^\)\(]*[a-z][^\)\(]*\)\.[a-z]',...
            'start','end','once');
          strucloc2 = strucloc2-1;
          [cellyloc1,cellyloc2] = regexpi(VarStr,'\{[^\}\{]*[a-z][^\}\{]*\}',...
            'start','end','once');
          if isempty(strucloc1) || (~isempty(cellyloc1) && strucloc1 > cellyloc1)
            strucloc1 = cellyloc1; strucloc2 = cellyloc2; svatype = 0;
          end
          while ~isempty(strucloc1)
            ADIGATOR.OPTIONS.PREALLOCATE = 1;
            % s(i).a = something; - have to replace this with
            % adigatorstrucasgn = something;
            % s = adigatorCellStrucAsgn(s,i,adigatorstrucasgn,'s(i).a','struct');
            fprintf(Tfid,[indent,'%% ',VarStr,'\n']);
            svacount = svacount+1;
            if svatype
              %svendloc = regexp(VarStr,'\([^\)]*[a-z][^\)]*\)\.\w*','end','once');
              svendloc = regexp(VarStr(strucloc2+1:end),'\W','once','start')+strucloc2-1;
              if isempty(svendloc)
                svendloc = length(VarStr);
              end
              StrucAsgn(svacount).i    = VarStr(strucloc1+1:strucloc2-2);
            else
              svendloc = strucloc2;
              StrucAsgn(svacount).i    = VarStr(strucloc1+1:strucloc2-1);
            end
            if DerNumber == 1
              StrucAsgn(svacount).new  = sprintf('adigatorStrucAsgn%1.0f',svacount);
            else
              StrucAsgn(svacount).new  = strtrim(FunStri(EqualLoc+1:end-1));
            end
            StrucAsgn(svacount).old  = VarStr(1:svendloc);
            StrucAsgn(svacount).s    = VarStr(1:strucloc1-1);
            StrucAsgn(svacount).type = svatype;
            VarStr  = [StrucAsgn(svacount).new,VarStr(svendloc+1:end)];
            VarStrings{Vcount} = VarStr;
            strucloc1 = strfind(FunStri,StrucAsgn(svacount).old);
            FunStri = [FunStri(1:strucloc1-1) StrucAsgn(svacount).new,...
              FunStri(strucloc1+length(StrucAsgn(svacount).old):end)];
            svatype = 1; % 1 for structure, 0 for cell
            [strucloc1,strucloc2] = regexpi(VarStr,'\([^\)\(]*[a-z][^\)\(]*\)\.[a-z]',...
              'start','end','once');
            strucloc2 = strucloc2-1;
            [cellyloc1,cellyloc2] = regexpi(VarStr,'\{[^\}]*[a-z][^\}]*\}',...
              'start','end','once');
            if isempty(strucloc1) || (~isempty(cellyloc1) && strucloc1 > cellyloc1)
              strucloc1 = cellyloc1; strucloc2 = cellyloc2; svatype = 0;
            end
          end
          
          % Check For SubsAsgn on all of the VarStrings
          if strcmp(VarStr(end),')')
            SubsFlags(Vcount) = 1;
            SubsLocs = strfind(VarStr,'(');
            StrucLoc = strfind(VarStr,').');
            if ~isempty(StrucLoc)
              SubsLocs = SubsLocs(SubsLocs > StrucLoc(end));
            end
            VarStrings{Vcount} = VarStr(1:SubsLocs(1)-1);
          end
        end
        
        if DerNumber == 1 || svacount == ADIGATOR.SVACOUNT
          % Print Out the Actual Calculation
          fprintf(Tfid,[indent,FunStri,'\n']);
          
          % Get the Strings for the variables that we are going to feed into
          % the Variable Analyzer
          if NumOutVars == 1
            LHSout = VarStrings{1};
            RHSout = sprintf([VarStrings{1},',''',VarStrings{1},...
              ''',%1.0f'],SubsFlags(1));
          else
            LHSout = cell(1,NumOutVars); RHSout = cell(1,NumOutVars);
            LHSout{1} = ['[',VarStrings{1},','];
            RHSout{1} = sprintf([VarStrings{1},',''',VarStrings{1},...
              ''',%1.0f,'],SubsFlags(1));
            for Vcount = 2:NumOutVars-1
              LHSout{Vcount} = [VarStrings{Vcount},','];
              RHSout{Vcount} = sprintf([VarStrings{Vcount},',''',VarStrings{Vcount},...
                ''',%1.0f,'],SubsFlags(Vcount));
            end
            LHSout{NumOutVars} = [VarStrings{NumOutVars},']'];
            RHSout{NumOutVars} = sprintf([VarStrings{NumOutVars},',''',VarStrings{NumOutVars},...
              ''',%1.0f'],SubsFlags(NumOutVars));
            LHSout = cell2mat(LHSout); RHSout = cell2mat(RHSout);
          end
          % Print Out the call to the Variable Analyzer
          FunStri = FindDoinkers(FunStri);
          fprintf(Tfid,[indent,LHSout,' = ',VAstr,'(''',FunStri,''',',RHSout,');\n']);
        end
        for svacounti = svacount:-1:ADIGATOR.SVACOUNT+1
          % Was a Structure Assignment
          fprintf(Tfid,[indent,'try adigatorDummyStruc = ',StrucAsgn(svacounti).s,...
            '; catch adigatorDummyErr; ']);
          if StrucAsgn(svacounti).type
            fprintf(Tfid,'adigatorDummyStruc = struct([]); end\n');
          else
            fprintf(Tfid,'adigatorDummyStruc = {}; end\n');
          end
          fprintf(Tfid,[indent,StrucAsgn(svacounti).s,' = adigatorCellStrucAsgn(%1.0f,',...
            'adigatorDummyStruc,',StrucAsgn(svacounti).new,...
            ',''',FindDoinkers(StrucAsgn(svacounti).old),' = ',StrucAsgn(svacounti).new,';'',',...
            StrucAsgn(svacounti).i,');\n'],...
            svacounti);
        end
        ADIGATOR.SVACOUNT = svacount;
      elseif StrLength > 7 && strcmp(FunStri(1:8),'keyboard')
        % KEYBOARD
        ADIGATOR.OPTIONS.KEYBOARD = 1;
        fprintf(Tfid,'keyboard\n');
        fprintf(Tfid,[indent,VAstr,'(''keyboard'');\n']);
      elseif StrLength > 4 && strcmp(FunStri(1:5),'break')
        % BREAK
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        fprintf(Tfid,[indent,'adigatorBreakCont(''break'',%1.0d,%1.0d);\n'],...
          IfCount,BroCount);
      elseif StrLength > 7 && strcmp(FunStri(1:8),'continue')
        % CONTINUE
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        fprintf(Tfid,[indent,'adigatorBreakCont(''continue'',%1.0d,%1.0d);\n'],...
          IfCount,BroCount);
      elseif StrLength >  6 && strcmp(FunStri(1:6),'error(')
        % ERROR MESSAGE
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        FunStri  = FindDoinkers(FunStri);
        fprintf(Tfid,[indent,'adigatorError(%1.0d,%1.0d,''',FunStri,''');\n'],...
          IfCount,BroCount);
      elseif StrLength > 5 && strcmp(FunStri(1:6),'return')
        % RETURN
        fprintf(Tfid,[indent,VAstr,'(''return'');\n']);
      elseif StrLength > 6 && strcmp(FunStri(1:7),'global ')
        % GLOBAL ASSIGNMENT
        VarString = strtrim(FunStri(8:end)); 
        if strcmp(VarString(end),';'); VarString = VarString(1:end-1); end
        VarLocs = strfind(VarString,' ');
        if ~isempty(VarLocs)
          % Check for multiple spaces.
          for Lcount = 2:length(VarLocs)
            if VarLocs(Lcount) == VarLocs(Lcount-1)+1
              VarLocs(Lcount-1) = 0;
            end
          end
          VarLocs = nonzeros(VarLocs);
          InVarStrs = cell(1,length(VarLocs)+1);
          OutVarStrs = cell(1,length(VarLocs)+1);
          InVarStrs{1} = ['''',strtrim(VarString(1:VarLocs(1)-1)),''','];
          OutVarStrs{1} = [strtrim(VarString(1:VarLocs(1)-1)),','];
          for Lcount = 2:length(VarLocs)
            InVarStrs{Lcount} = ['''',strtrim(VarString(...
              VarLocs(Lcount-1)+1:VarLocs(Lcount)-1)),''','];
            OutVarStrs{Lcount} = [strtrim(VarString(...
              VarLocs(Lcount-1)+1:VarLocs(Lcount)-1)),','];
          end
          InVarStrs{end} = ['''',strtrim(VarString(VarLocs(end)+1:end)),''''];
          OutVarStrs{end} = strtrim(VarString(VarLocs(end)+1:end));
          InVarStrs = cell2mat(InVarStrs);
          OutVarStrs = ['[',cell2mat(OutVarStrs),']'];
        else
          InVarStrs = ['''',VarString,''''];
          OutVarStrs = VarString;
        end
        fprintf(Tfid,[indent,OutVarStrs,' = ',VAstr,'(''global'',',InVarStrs,');\n']);
      elseif StrLength > 4 && strcmp(FunStri(1:5),'load(')
        % CHECK FOR LOAD - DONT ALLOW IT.
        errlink = GenErrorLink(Ffid,MajorLineCount);
        error(['In order to load in variables, must assign them to a'...
          ' variable name. At ',FunStri,' ',errlink]);
      else
        % Dont know what this statement is.
        errlink = GenErrorLink(Ffid,MajorLineCount);
        error(['Cannot process statement: ',FunStrFull,' at ',errlink])
      end
    end
  end
  MinorLineStart = 1;
  FunStrFULL = fgets(Ffid);
  MajorLineCount = MajorLineCount+1;
end

return
end

function FunStri = CheckFunctionCall(FunStri,FunStrChecks,Tfid,indent,DerNumber,MajorLineCount)

[Start,End] = regexp(FunStri,FunStrChecks,'start','end');

FunLoc = ~cellfun(@isempty,Start);
if FunLoc(1)
  errlink = GenErrorLink(Ffid,MajorLineCount);
  error(['User program cannot call main user function from within its methods at:' errlink])
elseif sum(FunLoc) > 1
  errlink = GenErrorLink(Ffid,MajorLineCount);
  error(['Currently can only handle one call to a user defined function per line at line: ',FunStri,' at: ',errlink])
elseif any(FunLoc)
  FunLocNum = 1:length(FunStrChecks);
  FunLoc = FunLocNum(FunLoc);
  Start = strfind(FunStri,FunStrChecks{FunLoc}(3:end));
  End   = End{FunLoc};
  if length(Start) > 1
    errlink = GenErrorLink(Ffid,MajorLineCount);
    error(['Currently can only handle one call to a user defined function per line at line: ',FunStri,' at: ',errlink])
  end
  fprintf(Tfid,[indent,'%% Call to User Function ',...
    FunStrChecks{FunLoc}(3:end-1),' --- (FunID %1.0d)\n'],FunLoc);
  % -------------------- Work on Inputs to Function --------------------- %
  InputStart = End+1;
  OpenLocs  = strfind(FunStri(InputStart:end),'(')+InputStart;
  CloseLocs = strfind(FunStri(InputStart:end),')')+InputStart;
  NumOpen   = length(OpenLocs);

  if ~isempty(OpenLocs)
      InputEnd = CloseLocs(find(CloseLocs(1:NumOpen) < OpenLocs,1,'first'))-2;
      if isempty(InputEnd)
        InputEnd = CloseLocs(NumOpen+1)-2;
      end
  else
    InputEnd = CloseLocs(1)-2;
  end
  InputStr = FunStri(InputStart:InputEnd);
  InputStrs = SeperateOutputStrings(InputStr(~isspace(InputStr)),1);
  NumInputs = length(InputStrs);
  InVarStrs = cell(1,NumInputs);
  for Icount = 1:NumInputs
    VarStr  = sprintf('cadainput%1.0d',Icount);
    InVarStrs{Icount} = [VarStr,';'];
    if DerNumber == 1 
      AsgnStr = [VarStr,' = ',InputStrs{Icount},';'];
      fprintf(Tfid,[indent,AsgnStr,'\n']);
      AsgnStr = FindDoinkers(AsgnStr);
      fprintf(Tfid,[indent,VarStr,' = adigatorVarAnalyzer(''',AsgnStr,''',',...
        VarStr,',''',VarStr,''',0);\n']);
    end
  end
  InVarStrs = cell2mat(InVarStrs);
  fprintf(Tfid,[indent,'adigatorInputs = {',InVarStrs(1:end-1),'};\n']);
  End = InputEnd+1;
  
  % -------------------- Print Call to Function ------------------------- %
  fprintf(Tfid,[indent,'[adigatorFunInfo, adigatorOutputs] = ',...
    'adigatortempfunc%1.0d(adigatorFunInfo,adigatorInputs);\n'],FunLoc);
  
  % ----------------------- Work On Outputs ----------------------------- %
  if ~isempty(regexp(FunStri(1:Start-1),'=\s*\>','once')) && strcmp(FunStri(1),'[') &&...
      (length(FunStri) == End+1 || ~isempty(regexp(FunStri(End+1:end),'<\\s*;','once')))
    % --------------------- Multiple Outputs ---------------------------- %
    EqLoc = strfind(FunStri,'=');
    OutputStart = 2;
    OutputEnd   = strfind(FunStri(1:EqLoc-1),']')-1;
    OutputStrs  = SeperateOutputStrings(FunStri(OutputStart:OutputEnd),0);
    NumOutput   = length(OutputStrs);
    for Ocount = 1:NumOutput
      OutVarStr = sprintf('cadaoutput%1.0d',Ocount);
      fprintf(Tfid,[indent,OutVarStr,' = adigatorOutputs{%1.0d};\n'],Ocount);
      if DerNumber == 1 && ~strcmp(OutputStrs{Ocount},'~')
        AsgnStr   = [OutputStrs{Ocount},' = ',OutVarStr,';'];
        fprintf(Tfid,[indent,AsgnStr,'\n']);
        SubLoc = strfind(OutputStrs{Ocount},'(');
        if ~isempty(SubLoc)
          SubsStr = '1';
          VarStr  = OutputStrs{Ocount}(1:SubLoc(1)-1);
        else
          SubsStr = '0';
          VarStr  = OutputStrs{Ocount};
        end
        fprintf(Tfid,[indent,VarStr,' = adigatorVarAnalyzer(''',AsgnStr,''',',...
          VarStr,',''',VarStr,''',',SubsStr,');\n']);
      end
    end
    FunStri = [];
  else
    % ----------------------- Single Output ----------------------------- %
    OutVarStr = 'cadaoutput1';
    fprintf(Tfid,[indent,OutVarStr,' = adigatorOutputs{1};\n']);
    if DerNumber == 1
      FunStri = [FunStri(1:Start-1),OutVarStr,FunStri(End+1:end)];
    else
      FunStri = [];
    end
  end
  
end
end

function FunStri = FindDoinkers(FunStri)

DoinkLocs = strfind(FunStri,'''');
if ~isempty(DoinkLocs)
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'''',FunStri(Dloc:end)];
  end
end
return
end

function FunStri = FindSlashes(FunStri,SlashLocs)

if ~isempty(SlashLocs)
  for Dcount = 1:length(SlashLocs)
    Dloc = SlashLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'\',FunStri(Dloc:end)];
  end
end
return
end

function FunStri = FindComments(FunStri)

DoinkLocs = strfind(FunStri,'%');
if ~isempty(DoinkLocs)
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'%',FunStri(Dloc:end)];
  end
end
return
end

function VarStrings = SeperateOutputStrings(VarStr,IOflag)
% Just seperates output variable strings
SpaceLocs = ones(1,length(VarStr));
SpaceLocs(~isspace(VarStr)) = 0;
CommaLocs = strfind(VarStr,',');
SpaceLocs(CommaLocs) = -1;
CharLocs  = 1:length(VarStr);

% Have a vector of locations where a zero or comma is.
ParenLoc1 = strfind(VarStr,'(');
if ~isempty(ParenLoc1)
  % Remove any entries of my zero/comma vector that are in
  % between parenthesis
  ParenLoc2 = strfind(VarStr,')');
  for Pcount = 1:length(ParenLoc1)
    SpaceLocs(CharLocs > ParenLoc1(Pcount) & CharLocs<ParenLoc2(Pcount)) = 0;
  end
end
CurlyLoc1 = strfind(VarStr,'{');
if ~isempty(CurlyLoc1)
  % Remove any entries of my zero/comma vector that are in
  % between curlies
  CurlyLoc2 = strfind(VarStr,'}');
  for Ccount = 1:length(CurlyLoc1)
    SpaceLocs(CharLocs > CurlyLoc1(Ccount) & CharLocs<CurlyLoc2(Ccount)) = 0;
  end
end
if IOflag
  SquareLoc1 = strfind(VarStr,'[');
  if ~isempty(SquareLoc1)
    % Remove any entries of my zero/comma vector that are in
    % between square brackets
    SquareLoc2 = strfind(VarStr,']');
    for Ccount = 1:length(SquareLoc1)
      SpaceLocs(SpaceLocs > SquareLoc1(Ccount) & SpaceLocs<SquareLoc2(Ccount)) = 0;
    end
  end
end
for Scount = 1:length(SpaceLocs)-1
  if SpaceLocs(Scount) && SpaceLocs(Scount+1)
    if SpaceLocs(Scount) < 1
      SpaceLocs(Scount+1) = 0;
      for S2count = Scount+2:length(SpaceLocs)-1
        if SpaceLocs(S2count) > 1
          SpaceLocs(S2count) = 1;
        else
          break
        end
      end
    else
      SpaceLocs(Scount) = 0;
    end
  end
end
SpaceLocs(end) = 0;
SpaceLocs1 = 1:length(VarStr);
SpaceLocs = SpaceLocs1(logical(SpaceLocs));

NumOutVars = length(SpaceLocs)+1;
VarLocs    = [0,SpaceLocs,length(VarStr)+1];
VarStrings = cell(NumOutVars,1);

for Vcount = 1:NumOutVars
  VarStrings{Vcount} = strtrim(VarStr(VarLocs(Vcount)+1:VarLocs(Vcount+1)-1));
end

end

function Statement = getIfForStatement(Ffid,Location)
fseek(Ffid,Location(4),-1);

StrFull = strtrim(fgets(Ffid));
while strcmp(StrFull(end-2:end),'...')
  StrFull = [StrFull(1:end-3),strtrim(fgets(Ffid))];
end

[MultStrs, ~] = adigatorSeperateFunLines(StrFull);
Statement     = MultStrs{Location(2)};
sLength       = length(Statement);
%Statement = FindDoinkers(Statement);
if sLength > 2 && strcmp(Statement(1:2),'if')
  Statement = strtrim(Statement(3:end));
elseif sLength > 6 && strcmp(Statement(1:6),'elseif')
  Statement = strtrim(Statement(7:end));
elseif sLength > 3 && strcmp(Statement(1:3),'for')
  Statement = strtrim(Statement(4:end));
else
  errlink = GenErrorLink(Ffid,Location(1));
  error(['??? Unsure of what this statement is:',Statement,' at: ',errlink])
end
Statement = regexprep(Statement,'&&','&');
Statement = regexprep(Statement,'\|\|','|');
if strcmp(Statement(end),';')
  Statement = strtrim(Statement(1:end-1));
end
if strcmp(Statement(end),',')
  Statement = strtrim(Statement(1:end-1));
end
end

function errlink = GenErrorLink(Ffid,MajorLineCount)
filename = fopen(Ffid);
filenameloc = strfind(filename,filesep);
filenameDisp = filename(filenameloc(end)+1:end);
errlink = sprintf(['User''s Function: <a href="matlab: opentoline(',filename,',%1.0d)">',filenameDisp,' line %1.0d</a>'],MajorLineCount,MajorLineCount);
end