function xOut = cadaOverMap(x,varargin)
%function varargout = cadaOverMap(x,varargin)
% this function Drives the Overmapping of Variables in For loops. If in the
% pre-printing run, it creates the OverMap by calling UnionVars. If in the
% printing run, it checks that the variable is properly overmapped. This is
% in the Overloaded folder because x is always overloaded.
% ------------------------ Input Information ---------------------------- %
% x:        variable to be overmapped in the FOR loop
% varargin: if there is a direct assignment, then this contains the .id
%           field of the assigned variable
% ------------------------ Output Information --------------------------- %
% xOut: variable to return to evaluating workspace
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR ADIGATORVARIABLESTORAGE
xOut = x;
if ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  xID  = x.id;
  if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF)
    % ------------------------ IF Overmap ------------------------------- %
    iOverLoc = ADIGATOR.VARINFO.OVERMAP.IF(xID,1);
    iOverForCount = ADIGATOR.VARINFO.OVERMAP.IF(xID,2);
    if iOverForCount
      % is a loop embedded within the if statement to which this overmap
      % belongs, need to check and see if this is the final iteration of
      % the loop.
      OverFlag = CheckForLoopIterations(iOverForCount,xID);
      if ~OverFlag
        iOverLoc = 0;
      end
    end
    if iOverLoc
      xOver = ADIGATORVARIABLESTORAGE.OVERMAP{iOverLoc};
      if ~isa(xOver,'cada') && ~ADIGATOR.EMPTYFLAG
        xOver = x;
      elseif ~ADIGATOR.EMPTYFLAG
        xOver = cadaUnionVars(xOver,x);
      end
      ADIGATORVARIABLESTORAGE.OVERMAP{iOverLoc} = xOver;
    end
    % ------------------------ RETURN ----------------------------------- %
    ReturnLoc = ADIGATOR.VARINFO.RETURN(xID);

    if ReturnLoc == 1 && iOverLoc
      % Return the OverMap
      xOut = xOver;
    else
      % Just Give back what came in.
      xOut = x;
    end
    % --------------------------- IF SAVE ------------------------------- %
    SaveLoc1 = ADIGATOR.VARINFO.SAVE.IF(xID,1);
    if SaveLoc1
      if iOverLoc
        ADIGATORVARIABLESTORAGE.SAVE{SaveLoc1} = xOver;
      else
        ADIGATORVARIABLESTORAGE.SAVE{SaveLoc1} = x;
      end
    end
  end
  
  if ~isempty(ADIGATOR.VARINFO.OVERMAP.FOR)
    SaveLoc = nonzeros(ADIGATOR.VARINFO.SAVE.FOR(xID,:));
    if ~isempty(SaveLoc)
      % Need to save the Variable for a FOR loop.
      if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF) &&...
          iOverLoc
        if isa(xOver,'cada')
          % Save the IF overmap
          ADIGATORVARIABLESTORAGE.SAVE{SaveLoc(1)} = xOver;
          if length(SaveLoc) == 2
            ADIGATORVARIABLESTORAGE.SAVE{SaveLoc(2)} = xOver;
          end
        end
      else
        % Save the regular var
        if isa(x,'cada')
          ADIGATORVARIABLESTORAGE.SAVE{SaveLoc(1)} = x;
          if length(SaveLoc) == 2
            ADIGATORVARIABLESTORAGE.SAVE{SaveLoc(2)} = x;
          end
        end
      end
    end
  end

  if ADIGATOR.FORINFO.FLAG
    % ------------------------ FOR LOOP OVERMAP ------------------------- %
    % We are in a FOR loop - Do the OverMapping (Note: this is separate
    % from the IF statement OverMapping.)
    if nargin == 2
      % Direct Assignment, we need to make sure that what is to the left of
      % the equal sign has the same overmapping as what is to the right of
      % the equals sign.
      DAid = varargin{1};
      if ~isempty(DAid)
        LHSloc = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
        RHSloc = ADIGATOR.VARINFO.OVERMAP.FOR(DAid,1);
        if LHSloc ~= RHSloc && RHSloc && LHSloc
          % If the RHS gets overmapped, then they share, if not (the RHS is
          % before the OverMap), then dont change anything.
          ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == LHSloc) = RHSloc;
          oldOverVar = ADIGATORVARIABLESTORAGE.OVERMAP{LHSloc};
          newOverVar = ADIGATORVARIABLESTORAGE.OVERMAP{RHSloc};
          if ~isempty(oldOverVar) && ~isempty(newOverVar)
            ADIGATORVARIABLESTORAGE.OVERMAP{RHSloc} = cadaUnionVars(oldOverVar,newOverVar);
          elseif ~isempty(oldOverVar)
            ADIGATORVARIABLESTORAGE.OVERMAP{RHSloc} = oldOverVar;
          end
%          elseif ~LHSloc && RHSloc
%            ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == RHSloc) = 0;
        end
      end
    end
    fOverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(xID);
    if fOverLoc
      xOver = ADIGATORVARIABLESTORAGE.OVERMAP{fOverLoc};
      if ~isa(xOver,'cada')
        xOver = x;
      else
        xOver = cadaUnionVars(xOver,x);
      end
      ADIGATORVARIABLESTORAGE.OVERMAP{fOverLoc} = xOver;
    end
    
    % -------------------- BREAK/CONTINUE OVERMAP ----------------------- %
    if ~isempty(ADIGATOR.VARINFO.OVERMAP.CONT)
      cOverLoc = ADIGATOR.VARINFO.OVERMAP.CONT(xID);
    else
      cOverLoc = 0;
    end
    if ~isempty(ADIGATOR.VARINFO.OVERMAP.BREAK)
      bOverLoc = ADIGATOR.VARINFO.OVERMAP.BREAK(xID);
    else
      bOverLoc = 0;
    end

    if (cOverLoc || bOverLoc) && ~ADIGATOR.EMPTYFLAG
      % Check to see if we are on last iteration of a loop or not
      InnerForLoc = ADIGATOR.FORINFO.INNERLOC;
      [ForLength, ForIter] = getForLength(InnerForLoc);
      if iOverLoc
        inVar = ADIGATORVARIABLESTORAGE.OVERMAP{iOverLoc};
      else
        inVar = x;
      end
      if ForIter < ForLength && cOverLoc
        % Do CONTINUE Overmap
        cOverVar = ADIGATORVARIABLESTORAGE.OVERMAP{cOverLoc};
        if isempty(cOverVar)
          ADIGATORVARIABLESTORAGE.OVERMAP{cOverLoc} = inVar;
        else
          ADIGATORVARIABLESTORAGE.OVERMAP{cOverLoc} = cadaUnionVars(cOverVar,inVar);
        end
      end
      if bOverLoc > 0 || (bOverLoc < 0 && ForIter == ForLength)
        % Do BREAK Overmap
        bOverLoc = abs(bOverLoc);
        bOverVar = ADIGATORVARIABLESTORAGE.OVERMAP{bOverLoc};
        if isempty(bOverVar)
          ADIGATORVARIABLESTORAGE.OVERMAP{bOverLoc} = inVar;
        else
          ADIGATORVARIABLESTORAGE.OVERMAP{bOverLoc} = cadaUnionVars(bOverVar,inVar);
        end
      end
    end
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  xID = x.id;
  if ADIGATOR.FORINFO.FLAG
    % ----------------------- In a FOR Loop ----------------------------- % 
    OverLoc  = ADIGATOR.VARINFO.OVERMAP.FOR(xID,1);
    SubsFlag = ADIGATOR.VARINFO.NAMELOCS(xID,3);
    if OverLoc
      xOver    = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
      if SubsFlag
        % Print the ReMap
        xOver = cadaPrintReMap(x,xOver,xID);
      else
        xOver = x;
      end
      
      if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF)
        % Check to see if we need to save for outside of conditional set
        IfSaveLoc = ADIGATOR.VARINFO.SAVE.IF(xID,1);
        if IfSaveLoc
          ADIGATORVARIABLESTORAGE.SAVE{IfSaveLoc(1)} = xOver;
        end
      end
      
      xOut = xOver;
      % Do Naming of the Variable
      funcname = cadafuncname(xID);
      xOut.func.name = funcname;
      for Vcount = 1:ADIGATOR.NVAROFDIFF
        if ~isempty(xOut.deriv(Vcount).nzlocs)
          xOut.deriv(Vcount).name = cadadername(funcname,Vcount,xID);
        end
      end
    end
  else
    % ---------------------- NOT in FOR Loop ---------------------------- %
    if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF)
      % Check to See if we Need to ReMap for INSIDE of conditional set
      OverLoc = ADIGATOR.VARINFO.OVERMAP.IF(xID,1);
      if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.FORINFO.OUTERLOC
        % Unrolling a loop which we are within, need to check to see if the
        % loop is contained within the conditional set 
      end
      
      % --------------------- IF REMAPPING ------------------------------ %
      if OverLoc
        xOver = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
        if isempty(xOver) && ~ADIGATOR.EMPTYFLAG
          xOver = x;
        elseif ~ADIGATOR.EMPTYFLAG && ~ADIGATOR.VARINFO.SAVE.IF(xID,2)
          xOver = cadaPrintReMap(x,xOver,xID);
        end
      end
      
      % ------------------------- RETURN -------------------------------- %
      ReturnLoc = ADIGATOR.VARINFO.RETURN(xID);
      if ReturnLoc
        % Return the OverMap
        xOut     = xOver;
        % Do Naming of the Variable
        funcname       = cadafuncname(xID);
        xOut.func.name = funcname;
        for Vcount = 1:ADIGATOR.NVAROFDIFF
          if ~isempty(xOut.deriv(Vcount).nzlocs)
            xOut.deriv(Vcount).name = cadadername(funcname,Vcount,xID);
          end
        end
      else
        % Just Give back what came in.
        xOut = x;
      end
      
      % --------------------------- IF SAVE ----------------------------- %
      SaveLoc1 = ADIGATOR.VARINFO.SAVE.IF(xID,1);
      if SaveLoc1
        if OverLoc
          ADIGATORVARIABLESTORAGE.SAVE{SaveLoc1} = xOver;
        else
          ADIGATORVARIABLESTORAGE.SAVE{SaveLoc1} = x;
        end
      end
      SaveLoc2 = ADIGATOR.VARINFO.SAVE.IF(xID,2);
      if SaveLoc2
        ADIGATORVARIABLESTORAGE.SAVE{SaveLoc2} = x;
      end
    end
    
    % Check to see if we need to Save for OUTSIDE of FOR loop
    if ~isempty(ADIGATOR.VARINFO.OVERMAP.FOR)
      ForSaveLoc = ADIGATOR.VARINFO.SAVE.FOR(xID,1);
      if ForSaveLoc
        % Need to save the Variable for a FOR loop.
        if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF) && OverLoc
          % Save the IF overmap
          ADIGATORVARIABLESTORAGE.SAVE{ForSaveLoc} = xOver;
        else
          % Save the regular var
          ADIGATORVARIABLESTORAGE.SAVE{ForSaveLoc} = x;
        end
      end
    end
  end
end

end

function flag = CheckForLoopIterations(ForLoc,xID)
global ADIGATORFORDATA
FORDATA = ADIGATORFORDATA;

WhileLoopFlag = 1;
while WhileLoopFlag
  ForLength = FORDATA(ForLoc).FOR(1).LENGTHS(end);
  ForIter   = FORDATA(ForLoc).COUNT.ITERATION;
  if ForLength ~= ForIter
    flag = 0;
    return
  end
  Children    = FORDATA(ForLoc).CHILDLOCS;
  if isempty(Children)
    flag = 1;
    return
  end
  ChildCounts = zeros(length(Children),2);
  for Ccount = 1:length(Children)
    Cloc = Children(Ccount);
    Cstart = FORDATA(Cloc).START;
    Cend   = FORDATA(Cloc).END;
    ChildCounts(Ccount,:) = [Cstart Cend];
  end
  ForLoc = Children(ChildCounts(:,1) <= xID & ChildCounts(:,2) >= xID);
  if isempty(ForLoc)
    flag = 1;
    break
  end
end
end

function [ForLength, ForIter] = getForLength(ForCount)
global ADIGATORFORDATA

ForLength   = ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(end);
ForIter     = ADIGATORFORDATA(ForCount).COUNT.ITERATION;

end