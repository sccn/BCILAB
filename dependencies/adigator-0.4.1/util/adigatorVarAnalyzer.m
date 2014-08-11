function varargout = adigatorVarAnalyzer(FunString,varargin)
%function varargout = adigatorVarAnalyzer(FunString,varargin)
%(non-overloaded)
% This module is the NON-Overloaded Version of adigatorVarAnalyzer. This is
% called from the temporary functions after a line of user code has
% been evaluated in order to analyze the outputs. As this is the non-overloaded
% version, it will only be called if all inputs are non-overloaded.
% -----------------------Input Information------------------------------- %
% FunString - the actual User's line of code which has just been evaluated
% varargin  - Information on all of the outputs from the user line of code
% which has just been evaluated. For each output, the inputs to
% adigatorVarAnalyzer are:
%         1. actual output
%         2. string of the name of the output (as defined by user)
%         3. flag stating whether output was subsasgn'd or not
% -----------------------Output Information------------------------------ %
% varargout - overloaded outputs, may or may not have some properties
% changed (may also be cells/structures of overloaded objects)
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
% ----------------------------- Parse Inputs ---------------------------- %
NUMvars   = nargout;
varargout = cell(NUMvars,1);
if strcmp(FunString,'global')
  % Global variables have a special case.
  VarStrings = varargin;
  Variables  = cadaGetGlobalVars(VarStrings);
  SubsFlags  = zeros(NUMvars,1);
else
  Variables  = cell(NUMvars,1);
  VarStrings = cell(NUMvars,1);
  SubsFlags  = zeros(NUMvars,1);
  for Vcount = 1:NUMvars
    Variables{Vcount}  = varargin{1+(Vcount-1)*3};
    VarStrings{Vcount} = varargin{2+(Vcount-1)*3};
    SubsFlags(Vcount)  = varargin{3+(Vcount-1)*3};
  end
end
if ADIGATOR.OPTIONS.PREALLOCATE
  varargout = Variables;
  return
end
PreOpCount = ADIGATOR.PREOPCOUNT;

if ADIGATOR.FILE.CALLFLAG
  % Function was called
  PreOpCount = adigatorFunctionCalled(Variables,VarStrings,SubsFlags);
  if ~PreOpCount
    % adigatorFunctionCalled took care of everything.
    ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
    ADIGATOR.FILE.CALLFLAG = 0;
    varargout = Variables;
    return
  end
  % Only way this wont fire as a return is that operations were done after
  % the function call. Ex. y = myfunc(x)*a;
end
if ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  if NUMvars
    % Is some sort of variable assignment.
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        % need to give this guy his own Operation
        % Count, and turn it into a symbolic
        x = adigatorMakeNumeric(x);
        ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1)...
          = ADIGATOR.VARINFO.COUNT;
        % Set his VarName
        adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,VarStrings{Vcount},SubsFlags(Vcount));
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      elseif isstruct(x) || iscell(x)
        % Is a cell or structure - Need to give a new operation count to
        % each numeric/overloaded object embedded in the cell/structure.
        % Numeric cells/fields need to be turned into Overloaded Objects.
        % Use adigatorFindCellStruc function (this case can occur in either
        % the overloaded version of cadaVarAnalyzer, or this version)
        if PreOpCount < ADIGATOR.VARINFO.COUNT && Vcount == 1
          xID = ADIGATOR.VARINFO.COUNT;
          ADIGATOR.VARINFO.NAMELOCS(PreOpCount:xID-1,2)=1:xID-PreOpCount;
        end
        if ADIGATOR.SUBSINDEXFLAG
          error(['??? Cannot Process Statement: ',FunString,...
            ' if using overloaded object as a reference of a structure/cell ,'...
            'result must be numeric or overloaded object']);
        elseif NUMvars > 1 && ~strcmp(FunString,'global')
          error(['??? Cannot Process Statement: ',FunString,...
            ' Multiple Outputs with structure/cell',...
            ' from a non-user function.']);
        end
        x = adigatorFindCellStruc(x,VarStrings{Vcount});
      else
%         error(['ADiGator only allows arrays, cells, and structures.',...
%           ' Including those objects embedded within cells/structures,',...
%           ' and global variables.']);
      end
      varargout{Vcount} = x;
    end
    if strcmp(FunString,'global')
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.PREOPCOUNT...
        :ADIGATOR.VARINFO.COUNT-1,3) = -Inf;
    end
  elseif strcmp(FunString,'keyboard')
    
  elseif strcmp(FunString,'return')
    
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if NUMvars
    % Is some sort of variable assignment.
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        % Variable is Numeric
        x = adigatorMakeNumeric(x);
        % OverMapping
        x = cadaOverMap(x);
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      elseif isstruct(x) || iscell(x)
        % Cell or Structure
        x = adigatorFindCellStruc(x,VarStrings{Vcount});
      end
      varargout{Vcount} = x;
    end
  elseif strcmp(FunString,'keyboard')
    
  elseif strcmp(FunString,'return')
    
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  if strcmp(FunString,'global')
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        x = adigatorMakeNumeric(x);
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      else
        x = adigatorFindCellStruc(x,VarStrings{Vcount});
      end
      VarStrings{Vcount} = [VarStrings{Vcount},' '];
      varargout{Vcount} = x;
    end
    fprintf(fid,[indent,'global ',cell2mat(VarStrings),'\n']);
  elseif ADIGATOR.FILE.CALLFLAG
    % A function was called - Everything will be worked out in
    % FunctionCatcher file.
    varargout = Variables;
  elseif NUMvars == 1
    % Single Variable Assignment - Either Numeric or Structure/Cell
    x = Variables{1};
    if isnumeric(x)
      % Need to give this its own Operation count and turn it into a
      % Symbolic.
      x = adigatorMakeNumeric(x);
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      % Can just print out whatever is to the right of the =
      if ~ADIGATOR.EMPTYFLAG
        EqLoc = strfind(FunString,'=');
        fprintf(fid,[indent,x.func.name,' = ',FunString(EqLoc+1:end),'\n']);
      end
      x = cadaOverMap(x);
    elseif (iscell(x) || isstruct(x))
      % Cell or Structure - Just print out the FunString and then call
      % adigatorFindCellStruc to change all of the variable names.
      if ~ADIGATOR.EMPTYFLAG
        FunString = FindComments(FunString);
        fprintf(fid,[indent,FunString,'\n']);
      end
      x = adigatorFindCellStruc(x,VarStrings{1});
    end
    varargout{Vcount} = x;
  elseif NUMvars > 1
    % Need to build what is to the left of the equals sign before we print
    % out anything. After we print it out, then we can call the OverMap In
    % this case, we should never have cells/structures because we didnt
    % allow in the trace run.
    LHSstrings = cell(1,NUMvars);
    for Vcount = 1:NUMvars
      % x is numeric.
      x = adigatorMakeNumeric(Variables{Vcount});
      if ~ADIGATOR.EMPTYFLAG
        LHSstrings{Vcount} = [x.func.name,','];
      end
      varargout{Vcount}  = x;
    end
    % Print out our LHS and the users RHS
    LHSstrings = cell2mat(LHSstrings);
    if ~ADIGATOR.EMPTYFLAG
      EQloc = strfind(FunString,'=');
      fprintf(fid,[indent,LHSstrings(1:end-1),' = ',...
        FunString(EQloc+1:end),'\n']);
    end
    % Make another run through the variables and check the overmapping.
    for Vcount = 1:NUMvar
      varargout{Vcount} = cadaOverMap(x);
    end
  elseif strcmp(FunString(1),'%')
    % User Comment.
    FunString = FindComments(FunString);
  elseif strcmp(FunString,'keyboard')
    % KEYBOARD
  elseif strcmp(FunString,'return')
    % RETURN STATEMENT
  elseif strfind(FunString,'global')
    
  end

  if ADIGATOR.DERNUMBER == 1 && ADIGATOR.OPTIONS.COMMENTS
    fprintf(fid,[indent,'%%User Line: ',FunString,'\n']);
  elseif length(FunString) > 1 && strcmp(FunString(1),'%')
    fprintf(fid,[indent,FunString,'\n']);
  elseif ADIGATOR.OPTIONS.COMMENTS
    fprintf(fid,[indent,'%% Deriv %1.0d Line: ',...
      FunString,'\n'],ADIGATOR.DERNUMBER-1);
  end
end

ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
ADIGATOR.FILE.CALLFLAG = 0;
ADIGATOR.SUBSINDEXFLAG = 0;
return
end

function FunStri = FindComments(FunStri)

DoinkLocs = strfind(FunStri,'%');
if ~isempty(DoinkLocs)
  NUMdoinks = length(DoinkLocs);
  FunStri(end+NUMdoinks) = ' ';
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'%',FunStri(Dloc:end-1)];
  end
end
return
end

function cadaVars = cadaGetGlobalVars(cadaGlobalStrs)

cadaVars = cell(size(cadaGlobalStrs));
for cadaVcount = 1:length(cadaGlobalStrs)
  eval(['global ',cadaGlobalStrs{cadaVcount},';']);
  cadaVars{cadaVcount} = eval(cadaGlobalStrs{cadaVcount});
end

end