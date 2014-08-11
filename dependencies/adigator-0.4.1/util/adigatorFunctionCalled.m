function flag = adigatorFunctionCalled(Variables,VarStrings,SubsFlags)
% function flag = adigatorFunctionCalled(Variables,VarStrings,SubsFlags)
% This function is called whenever we get out of a function call. This
% needs to look and check for anything happening AFTER the function call.
% Was there an operation performed after the function call? Was it just a
% subsasgn? Was it like y = x*myfunc(x)?
% ------------------------ Input Information ---------------------------- %
% Variables:  actual outputs from the function evaluation
% VarStrings: user defined string names of the outputs
% SubsFlags:  flags on which variables are subsasgned on the output
% ------------------------ Output Information --------------------------- %
% flag: flag telling adigatorVarAnalyzer if there were operations after the
% function call
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFUNCTIONINFO
FunID = ADIGATOR.FILE.CALLFLAG;

if ADIGATOR.FORINFO.FLAG == 1
  LastVar = ADIGATORFUNCTIONINFO(FunID).OUTPUT.EMPTYSTRUCVARS{end};
else
  CallCount = ADIGATORFUNCTIONINFO(FunID).ITERATION.CALLCOUNT;
  IterCount = ADIGATORFUNCTIONINFO(FunID).ITERATION.ITERID(CallCount);
  LastVar   = ADIGATORFUNCTIONINFO(FunID).OUTPUT.STRUCVARS{end,IterCount};
end
LastFuncOp = LastVar.id;


if LastFuncOp < ADIGATOR.VARINFO.COUNT-1 &&...
    sum(SubsFlags) < ADIGATOR.VARINFO.COUNT-1-LastFuncOp
  % The operations that took place were something else, so the output of
  % the function was not directly assigned.
  % Ex. y = x.*myfunc(x);
  % Set up some temporary structures to be the outputs of the function
  % call
  if ADIGATOR.FORINFO.FLAG < 2
    TempVarCount = 0;
    if ADIGATOR.FORINFO.FLAG == 1
      FunctionOutputs = ADIGATORFUNCTIONINFO(FunID)...
        .OUTPUT.EMPTYSTRUCVARS;
    else
      FunctionOutputs = ADIGATORFUNCTIONINFO(FunID)...
        .OUTPUT.STRUCVARS(:,IterCount);
    end
    for Ocount = 1:length(FunctionOutputs)
      TempVarCount = TempVarCount+1;
      TempVarName = sprintf('cada%1.0dtemps%1.0d',...
        ADIGATOR.DERNUMBER,TempVarCount);
      OpCount = FunctionOutputs{Ocount}.id;
      adigatorAssignImpVarNames(OpCount,TempVarName,0);
      ADIGATOR.VARINFO.NAMELOCS(OpCount,2) = 0;
    end
  end
  
  % Set flag to tell cadaVarAnalyzer to go about analyzing the other
  % outputs as it normally would.
  flag = LastFuncOp+1;
else
  % The Operations that took place were just SubsAsgns into the outputs.
  % Ex. [y(1), y(2)] = myfunc(x) OR no operations took place.
  % Ex. y = myfunc(x);
  if ADIGATOR.FORINFO.FLAG < 2
    TempVarCount = 0;
    FunctionNames     = ADIGATORFUNCTIONINFO(FunID).OUTPUT.NAMES;
    FunctionNamesFull = ADIGATORFUNCTIONINFO(FunID).OUTPUT.STRUCNAMES;
    for Vcount = 1:length(Variables)
      if SubsFlags(Vcount)
        % SubsAsgn case - need to get the operation count of the variable
        % going into SubsAsgn
        OutName = FunctionNames{Vcount};
        OutLoc = strcmp(OutName,FunctionNamesFull);
        if ADIGATOR.FORINFO.FLAG == 1
          OutVar = ADIGATORFUNCTIONINFO(FunID)...
            .OUTPUT.EMPTYVARS{OutLoc};
        else
          OutVar = ADIGATORFUNCTIONINFO(FunID)...
            .OUTPUT.VARS{OutLoc,IterCount};
        end
        OutOp = OutVar.id;
        % Assign the variable going into subsasgn as an important variable
        TempVarCount = TempVarCount+1;
        TempVarName = sprintf('cada%1.0dtemps%1.0d',...
          ADIGATOR.DERNUMBER,TempVarCount);
        adigatorAssignImpVarNames(OutOp,TempVarName,0);
        
        % Need to assign the variable name of the output operation too.
        OpCount = Variables{Vcount}.id;
        adigatorAssignImpVarNames(OpCount,VarStrings{Vcount},1);
      elseif isa(Variables{Vcount},'cada')
        % Just a normal overloaded object.
        OpCount = Variables{Vcount}.id;
        adigatorAssignImpVarNames(OpCount,VarStrings{Vcount},0);
      else
        % Cell or Structure
        % We already have the info of everything buried in it saved in
        % ADIGATORFUNCTIONINFO - just need to get it.
        OutName   = FunctionNames{Vcount};
        OutLength = length(OutName);
        for Scount = 1:length(FunctionNamesFull)
          SOutName = FunctionNamesFull{Scount};
          if length(SOutName) > OutLength &&...
              strcmp(OutName,SOutName(1:OutLength))
            if ADIGATOR.FORINFO.FLAG == 1
              OutVar = ADIGATORFUNCTIONINFO(FunID)...
                .OUTPUT.EMPTYSTRUCVARS{Scount};
            else
              OutVar = ADIGATORFUNCTIONINFO(FunID)...
                .OUTPUT.STRUCVARS{Scount,IterCount};
            end
            OpCount = OutVar.id;
            OutVarStr = [VarStrings{Vcount},SOutName(OutLength+1:end)];
            adigatorAssignImpVarNames(OpCount,OutVarStr,0);
            if ADIGATOR.DERNUMBER == 1 && ADIGATORFUNCTIONINFO(FunID).DERNUMBER > 1
              ADIGATOR.VARINFO.NAMELOCS(OpCount,2) = inf;
              ADIGATOR.VARINFO.NAMELOCS(OpCount,3) = FunID;
            end
          end
        end
      end
    end
  end
  flag = 0;
end