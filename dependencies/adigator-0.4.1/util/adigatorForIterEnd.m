function [outEvalStr outEvalVar] = adigatorForIterEnd(ForCount,ForIter)
% This transformation routine is placed at the end of each iteration of a
% FOR loop in the intermediate program.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFORDATA

if ADIGATOR.OPTIONS.PREALLOCATE
 % Pre-Allocation run for structures/cells
 outEvalVar = [];
 outEvalStr = [];
elseif ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  ADIGATORFORDATA(ForCount).END = ADIGATOR.VARINFO.COUNT-1;
  if ~ADIGATOR.OPTIONS.UNROLL
    % ------------------ Pre-Allocate ADIGATORFORDATA ---------------------- %
    NUMPref   = ADIGATORFORDATA(ForCount).COUNT.SUBSREF;
    % Number of SUBREFS in this loop
    NUMPasgn  = ADIGATORFORDATA(ForCount).COUNT.SUBSASGN;
    % Number of SUBSASGNs in this loop
    NUMPsp    = ADIGATORFORDATA(ForCount).COUNT.SPARSE;
    % Number of SPARSEs in this loop
    NUMPnz    = ADIGATORFORDATA(ForCount).COUNT.NONZEROS;
    % Number of NONZEROS in this loop
    NUMPhorz  = ADIGATORFORDATA(ForCount).COUNT.HORZCAT;
    % Number of HORZCATS in this loop
    NUMPvert  = ADIGATORFORDATA(ForCount).COUNT.VERTCAT;
    % Number of VERTCATS in this loop
    NUMPtran  = ADIGATORFORDATA(ForCount).COUNT.TRANSPOSE;
    % Number of TRANSPOSES in this loop
    NUMPrep   = ADIGATORFORDATA(ForCount).COUNT.REPMAT;
    % Number of REPMATs in this loop
    NUMPres   = ADIGATORFORDATA(ForCount).COUNT.RESHAPE;
    % Number of RESHAPEs in this loop
    NUMPsiz   = ADIGATORFORDATA(ForCount).COUNT.SIZE;
    % Number of SIZESs in this loop.
    NUMchild = size(ADIGATORFORDATA(ForCount).CHILDLOCS,1);
    % Number of Children in this loop
    NUMTref   = NUMPref;   NUMTasgn  = NUMPasgn;  NUMTfor  = 1;
    NUMTnz    = NUMPnz;    NUMThorz  = NUMPhorz;  NUMTvert = NUMPvert;
    NUMTsp    = NUMPsp;    NUMTrep   = NUMPrep;   NUMTres  = NUMPres;
    NUMTtran  = NUMPtran;  NUMTsiz   = NUMPsiz;
    for Ccount = 1:NUMchild
      % Get the total number of Special Function Calls
      ChildLoc  = ADIGATORFORDATA(ForCount).CHILDLOCS(Ccount,1);
      NUMTref   = length(ADIGATORFORDATA(ChildLoc).SUBSREF)   + NUMTref;
      NUMTasgn  = length(ADIGATORFORDATA(ChildLoc).SUBSASGN)  + NUMTasgn;
      NUMTsp    = length(ADIGATORFORDATA(ChildLoc).SPARSE)    + NUMTsp;
      NUMTnz    = length(ADIGATORFORDATA(ChildLoc).NONZEROS)  + NUMTnz;
      NUMTfor   = length(ADIGATORFORDATA(ChildLoc).FOR)       + NUMTfor;
      NUMThorz  = length(ADIGATORFORDATA(ChildLoc).HORZCAT)   + NUMThorz;
      NUMTvert  = length(ADIGATORFORDATA(ChildLoc).VERTCAT)   + NUMTvert;
      NUMTtran  = length(ADIGATORFORDATA(ChildLoc).TRANSPOSE) + NUMTtran;
      NUMTrep   = length(ADIGATORFORDATA(ChildLoc).REPMAT)    + NUMTrep;
      NUMTres   = length(ADIGATORFORDATA(ChildLoc).RESHAPE)   + NUMTres;
      NUMTsiz   = length(ADIGATORFORDATA(ChildLoc).SIZE)      + NUMTsiz;
    end
    
    if NUMTref
      % Pre-Allocation - SUBSREF
      ADIGATORFORDATA(ForCount).SUBSREF = ...
        struct('VARS',[],'FLAGS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).SUBSREF(NUMTref,1).VARS = [];
    end
    if NUMTasgn
      % Pre-Allocation - SUBSASGN
      ADIGATORFORDATA(ForCount).SUBSASGN = ...
        struct('VARS',[],'FLAGS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).SUBSASGN(NUMTasgn,1).VARS = [];
    end
    if NUMTnz
      % Pre-Allocation - NONZEROS
      ADIGATORFORDATA(ForCount).NONZEROS = ...
        struct('VARS',[],'FLAGS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).NONZEROS(NUMTnz,1).VARS = [];
    end
    if NUMTsp
      % Pre-Allocation - SPARSE
      ADIGATORFORDATA(ForCount).SPARSE = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).SPARSE(NUMTsp,1).VARS = [];
    end
    
    if NUMThorz
      % Pre-Allocation - HORZCAT
      ADIGATORFORDATA(ForCount).HORZCAT = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).HORZCAT(NUMThorz,1).VARS = [];
    end
    if NUMTvert
      % Pre-Allocation - VERTCAT
      ADIGATORFORDATA(ForCount).VERTCAT = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).VERTCAT(NUMTvert,1).VARS = [];
    end
    if NUMTtran
      % Pre-Allocation - TRANSPOSE
      ADIGATORFORDATA(ForCount).TRANSPOSE = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).TRANSPOSE(NUMTtran,1).VARS = [];
    end
    if NUMTrep
      % Pre-Allocation - REPMAT
      ADIGATORFORDATA(ForCount).REPMAT = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).REPMAT(NUMTrep,1).VARS = [];
    end
    if NUMTres
      % Pre-Allocation - RESHAPE
      ADIGATORFORDATA(ForCount).RESHAPE = ...
        struct('VARS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).RESHAPE(NUMTres,1).VARS = [];
    end
    if NUMTsiz
      % Pre-Allocation - SIZE
      ADIGATORFORDATA(ForCount).SIZE = ...
        struct('SIZES',[],'LOCS',[]);
      ADIGATORFORDATA(ForCount).SIZE(NUMTsiz,1).SIZES = [];
    end
    % Pre-Allocation - CHILD LENGTHS
    ADIGATORFORDATA(ForCount).FOR = ...
      struct('LENGTHS',[],'LOCS',[]);
    ADIGATORFORDATA(ForCount).FOR(NUMTfor,1).LENGTHS = [];
    % Pre-Allocation - OTHER
    NUMother = ADIGATORFORDATA(ForCount).COUNT.OTHER;
    if NUMother
      ADIGATORFORDATA(ForCount).OTHER = ...
        struct('DATA',[]);
      ADIGATORFORDATA(ForCount).OTHER(NUMother,1).DATA = [];
    end
    %------------- Organizational Functions in Parent Loop -----------------%
    % SUBSREF locations in this (parent) loop
    for Rcount = 1:NUMPref
      ADIGATORFORDATA(ForCount).SUBSREF(Rcount).LOCS   =   [ForCount;Rcount];
    end
    Rcount = NUMPref;
    % SUBSASGN locations in this (parent) loop
    for Scount = 1:NUMPasgn
      ADIGATORFORDATA(ForCount).SUBSASGN(Scount).LOCS  =   [ForCount;Scount];
    end
    Scount = NUMPasgn;
    % SPARSE locations in this (parent) loop
    for SPcount = 1:NUMPsp
      ADIGATORFORDATA(ForCount).SPARSE(SPcount).LOCS   =   [ForCount;SPcount];
    end
    SPcount = NUMPsp;
    % NONZEROS locations in this (parent) loop
    for Ncount = 1:NUMPnz
      ADIGATORFORDATA(ForCount).NONZEROS(Ncount).LOCS  =   [ForCount;Ncount];
    end
    Ncount = NUMPnz;
    % HORZCAT locations in this (parent) loop
    for Hcount = 1:NUMPhorz
      ADIGATORFORDATA(ForCount).HORZCAT(Hcount).LOCS   =   [ForCount;Hcount];
    end
    Hcount = NUMPhorz;
    % VERTCAT locations in this (parent) loop
    for Vcount = 1:NUMPvert
      ADIGATORFORDATA(ForCount).VERTCAT(Vcount).LOCS   =   [ForCount;Vcount];
    end
    Vcount = NUMPvert;
    % TRANSPOSE locations in this (parent) loop
    for Tcount = 1:NUMPtran
      ADIGATORFORDATA(ForCount).TRANSPOSE(Tcount).LOCS =   [ForCount;Tcount];
    end
    Tcount = NUMPtran;
    % REPMAT locations in this (parent) loop
    for REPcount = 1:NUMPrep
      ADIGATORFORDATA(ForCount).REPMAT(REPcount).LOCS  =   [ForCount;REPcount];
    end
    REPcount = NUMPrep;
    % RESHAPE locations in this (parent) loop
    for REScount = 1:NUMPres
      ADIGATORFORDATA(ForCount).RESHAPE(REScount).LOCS =   [ForCount;REScount];
    end
    REScount = NUMPres;
    % SIZE locations in this (parent) loop
    for SIZcount = 1:NUMPsiz
      ADIGATORFORDATA(ForCount).SIZE(SIZcount).LOCS    =   [ForCount;SIZcount];
    end
    SIZcount = NUMPsiz;
    % FORLENGTH location in this (parent) loop
    ADIGATORFORDATA(ForCount).FOR(1).LOCS = [ForCount;1];
    Fcount = 1;
    
    % --------- Organizational Function Locations in Child Loops -----------%
    for Ccount = 1:NUMchild
      ChildLoc = ADIGATORFORDATA(ForCount).CHILDLOCS(Ccount,1);
      CHILDDATA = ADIGATORFORDATA(ChildLoc);
      for R2count = 1:length(CHILDDATA.SUBSREF)
        Rcount = Rcount+1;
        ADIGATORFORDATA(ForCount).SUBSREF(Rcount).LOCS = ...
          [[ForCount;Rcount],CHILDDATA.SUBSREF(R2count).LOCS];
      end
      for S2count = 1:length(CHILDDATA.SUBSASGN)
        Scount = Scount+1;
        ADIGATORFORDATA(ForCount).SUBSASGN(Scount).LOCS = ...
          [[ForCount;Scount],CHILDDATA.SUBSASGN(S2count).LOCS];
      end
      for SP2count = 1:length(CHILDDATA.SPARSE)
        SPcount = SPcount+1;
        ADIGATORFORDATA(ForCount).SPARSE(SPcount).LOCS = ...
          [[ForCount;SPcount],CHILDDATA.SPARSE(SP2count).LOCS];
      end
      for N2count = 1:length(CHILDDATA.NONZEROS)
        Ncount = Ncount+1;
        ADIGATORFORDATA(ForCount).NONZEROS(Ncount).LOCS = ...
          [[ForCount;Ncount],CHILDDATA.NONZEROS(N2count).LOCS];
      end
      for H2count = 1:length(CHILDDATA.HORZCAT)
        Hcount = Hcount+1;
        ADIGATORFORDATA(ForCount).HORZCAT(Hcount).LOCS = ...
          [[ForCount;Hcount],CHILDDATA.HORZCAT(H2count).LOCS];
      end
      for V2count = 1:length(CHILDDATA.VERTCAT)
        Vcount = Vcount+1;
        ADIGATORFORDATA(ForCount).VERTCAT(Vcount).LOCS = ...
          [[ForCount;Vcount],CHILDDATA.VERTCAT(V2count).LOCS];
      end
      for T2count = 1:length(CHILDDATA.TRANSPOSE)
        Tcount = Tcount+1;
        ADIGATORFORDATA(ForCount).TRANSPOSE(Tcount).LOCS = ...
          [[ForCount;Tcount],CHILDDATA.TRANSPOSE(T2count).LOCS];
      end
      for REP2count = 1:length(CHILDDATA.REPMAT)
        REPcount = REPcount+1;
        ADIGATORFORDATA(ForCount).REPMAT(REPcount).LOCS = ...
          [[ForCount;REPcount],CHILDDATA.REPMAT(REP2count).LOCS];
      end
      for RES2count = 1:length(CHILDDATA.RESHAPE)
        REScount = REScount+1;
        ADIGATORFORDATA(ForCount).RESHAPE(REScount).LOCS = ...
          [[ForCount;REScount],CHILDDATA.RESHAPE(RES2count).LOCS];
      end
      for SIZ2count = 1:length(CHILDDATA.SIZE)
        SIZcount = SIZcount+1;
        ADIGATORFORDATA(ForCount).SIZE(SIZcount).LOCS = ...
          [[ForCount;SIZcount],CHILDDATA.SIZE(SIZ2count).LOCS];
      end
      for F2count = 1:length(CHILDDATA.FOR)
        Fcount = Fcount+1;
        ADIGATORFORDATA(ForCount).FOR(Fcount).LOCS = ...
          [[ForCount;Fcount],CHILDDATA.FOR(F2count).LOCS];
      end
    end
  end
  
  % --------------------- Set up Exit of FOR Loop ----------------------- %
  PARENTLOC = ADIGATORFORDATA(ForCount).PARENTLOC;
  if PARENTLOC
    % Still within a FOR loop
    ADIGATOR.FORINFO.INNERLOC = PARENTLOC;
  else
    % This is parent FOR loop
    ADIGATOR.FORINFO.INNERLOC = 0;
    ADIGATOR.FORINFO.OUTERLOC = 0;
    ADIGATOR.FORINFO.FLAG     = 0;
  end
  outEvalStr = [];
  outEvalVar = [];
elseif ADIGATOR.OPTIONS.UNROLL
  % --------------------------------------------------------------------- %
  %                           Unrolling the Loop                          %
  % --------------------------------------------------------------------- %
  ForLength = ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(end);
  if ForIter == ForLength
    % Last Iteration
    PARENTLOC = ADIGATORFORDATA(ForCount).PARENTLOC;
    if PARENTLOC
      % Still within a FOR loop
      ADIGATOR.FORINFO.INNERLOC = PARENTLOC;
      ADIGATORFORDATA(PARENTLOC).CHILDLOCS(end+1,1) = ForCount;
    else
      % This is parent FOR loop
      ADIGATOR.FORINFO.INNERLOC = 0;
      ADIGATOR.FORINFO.OUTERLOC = 0;
      ADIGATOR.FORINFO.FLAG     = 0;
    end
  end
  outEvalStr = [];
  outEvalVar = [];
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  
  % --------------------- Gather Children Data -------------------------- %
  if ~isempty(ADIGATORFORDATA(ForCount).CHILDLOCS)
    % ---------GATHER CHILDREN STUFF----------------- %
    % Gather Children SUBSREF Indices
    for Rcount = ADIGATORFORDATA(ForCount).COUNT.SUBSREF+1:...
        length(ADIGATORFORDATA(ForCount).SUBSREF)
      CollectChildData(ForCount,ForIter,Rcount,'SUBSREF','INDICES');
    end
    % Gather Children SUBSASGN Indices
    for Scount = ADIGATORFORDATA(ForCount).COUNT.SUBSASGN+1:...
        length(ADIGATORFORDATA(ForCount).SUBSASGN)
      CollectChildData(ForCount,ForIter,Scount,'SUBSASGN','INDICES');
    end
    % Gather Children SPARSE Indices
    for SPcount = ADIGATORFORDATA(ForCount).COUNT.SPARSE+1:...
        length(ADIGATORFORDATA(ForCount).SPARSE)
      CollectChildData(ForCount,ForIter,SPcount,'SPARSE','INDICES');
    end
    % Gather Children NONZEROS Indices
    for Ncount = ADIGATORFORDATA(ForCount).COUNT.NONZEROS+1:...
        length(ADIGATORFORDATA(ForCount).NONZEROS)
      CollectChildData(ForCount,ForIter,Ncount,'NONZEROS','INDICES');
    end
    % Gather Children HORZCAT Sizes
    for Hcount = ADIGATORFORDATA(ForCount).COUNT.HORZCAT+1:...
        length(ADIGATORFORDATA(ForCount).HORZCAT)
      CollectChildData(ForCount,ForIter,Hcount,'HORZCAT','SIZES');
    end
    % Gather Children VERTCAT Sizes
    for Vcount = ADIGATORFORDATA(ForCount).COUNT.VERTCAT+1:...
        length(ADIGATORFORDATA(ForCount).VERTCAT)
      CollectChildData(ForCount,ForIter,Vcount,'VERTCAT','SIZES');
    end
    % Gather Children TRANSPOSE Sizes
    for Tcount = ADIGATORFORDATA(ForCount).COUNT.TRANSPOSE+1:...
        length(ADIGATORFORDATA(ForCount).TRANSPOSE)
      CollectChildData(ForCount,ForIter,Tcount,'TRANSPOSE','SIZES');
    end
    % Gather Children REPMAT Sizes
    for REPcount = ADIGATORFORDATA(ForCount).COUNT.REPMAT+1:...
        length(ADIGATORFORDATA(ForCount).REPMAT)
      CollectChildData(ForCount,ForIter,REPcount,'REPMAT','SIZES');
    end
    % Gather Children RESHAPE Sizes
    for REScount = ADIGATORFORDATA(ForCount).COUNT.RESHAPE+1:...
        length(ADIGATORFORDATA(ForCount).RESHAPE)
      CollectChildData(ForCount,ForIter,REScount,'RESHAPE','SIZES');
    end
    % Gather Children SIZE Sizes
    for SIZcount = ADIGATORFORDATA(ForCount).COUNT.SIZE+1:...
        length(ADIGATORFORDATA(ForCount).SIZE)
      CollectChildData(ForCount,ForIter,SIZcount,'SIZE','SIZES');
    end
    % Gather Children FORLENGTHS
    for Fcount = 2:length(ADIGATORFORDATA(ForCount).FOR)
      CFloc     = ADIGATORFORDATA(ForCount).FOR(Fcount).LOCS(1,2);
      CFcount   = ADIGATORFORDATA(ForCount).FOR(Fcount).LOCS(2,2);
      if CFcount == 1
        % Since we actually store forlengths in the top slot across the
        % parents iterations, we only need the last one.
        CFlengths = ADIGATORFORDATA(CFloc).FOR(CFcount).LENGTHS(end);
      else
        CFlengths = ADIGATORFORDATA(CFloc).FOR(CFcount).LENGTHS;
      end
      NUMCF     = numel(CFlengths);
      
      if ForIter > 1
        ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS...
          (1:NUMCF,ForIter) = CFlengths(:);
      else
        ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS = CFlengths(:);
      end
    end
  end
  
  % --------------------- See if on Last Iter --------------------------- %
  ForLength = ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(end);
  if ForIter == ForLength
    % Last Iteration
    PARENTLOC = ADIGATORFORDATA(ForCount).PARENTLOC;
    if PARENTLOC
      % Still within a FOR loop
      ADIGATOR.FORINFO.INNERLOC = PARENTLOC;
      ADIGATORFORDATA(PARENTLOC).CHILDLOCS(end+1,1) = ForCount;
    else
      % This is parent FOR loop
      ADIGATOR.FORINFO.INNERLOC = 0;
      ADIGATOR.FORINFO.OUTERLOC = 0;
      ADIGATOR.FORINFO.FLAG     = 0;
    end
  end
  
  % ---------------------- Returning BREAK/CONTINUE Stuff --------------- %
  if ForIter < ForLength && ~isempty(ADIGATOR.VARINFO.OVERMAP.CONT)
    % We need to return the CONTINUE OVERMAPs associated with this loop to
    % the evaluating workspace
    Start = ADIGATORFORDATA(ForCount).START;
    End   = ADIGATORFORDATA(ForCount).END;
    LoopCounts = Start:End;
    ContOverLocs  = nonzeros(unique(ADIGATOR.VARINFO.OVERMAP.CONT(LoopCounts)));
    [outEvalStr, outEvalVar]  = getContBreakOvermaps(ContOverLocs,1,LoopCounts);
  elseif ForIter == ForLength && ~isempty(ADIGATOR.VARINFO.OVERMAP.BREAK)
    Start = ADIGATORFORDATA(ForCount).START;
    End   = ADIGATORFORDATA(ForCount).END;
    LoopCounts = Start:End;
    BreakOverLocs  = nonzeros(unique(abs(ADIGATOR.VARINFO.OVERMAP.BREAK(LoopCounts))));
    [outEvalStr, outEvalVar]  = getContBreakOvermaps(BreakOverLocs,0,LoopCounts);
  else
    outEvalStr = [];
    outEvalVar = [];
  end
elseif ADIGATOR.RUNFLAG == 2
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;
  if ParentLoc
    ADIGATOR.FORINFO.INNERLOC = ParentLoc;
    outEvalStr = [];
    outEvalVar = [];
  else
    ADIGATOR.FORINFO.INNERLOC = 0;
    ADIGATOR.FORINFO.OUTERLOC = 0;
    ADIGATOR.FORINFO.FLAG     = 0;
  end
  fid    = ADIGATOR.PRINT.FID;
  ADIGATOR.PRINT.INDENT(1:4) = [];
  indent = ADIGATOR.PRINT.INDENT;
  if ~ADIGATOR.EMPTYFLAG
    fprintf(fid,[indent,'end\n']);
  end
  if ~ParentLoc
    % ---------------- Parent Loop - Do ReMapping ----------------------- %
    VarCounts = ADIGATORFORDATA(ForCount).START:ADIGATORFORDATA(ForCount).END;
    [outEvalStr, outEvalVar] = DoRemapping(VarCounts,ForCount);
  end
end
end

function CollectChildData(ForCount,ForIter,SFcount,FunStr,RefStr)
global ADIGATORFORDATA
% Fuction is used to collect data from Children. SFcount (Special Function
% Count) is the count in the parent (ForCount). FunStr is the string
% corresponding to this special function (i.e. 'SUBSREF', 'SUBSASGN',etc.)
% RefStr is the string corresponding to where the data is being collected
% from (either 'INDICES' or 'SIZES')

CSFloc    = ADIGATORFORDATA(ForCount).(FunStr)(SFcount).LOCS(1,2);
CSFcount  = ADIGATORFORDATA(ForCount).(FunStr)(SFcount).LOCS(2,2);
CSFdata   = ADIGATORFORDATA(CSFloc).  (FunStr)(CSFcount).(RefStr)(:);
OldData   = ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr);
NUMCSF    = length(CSFdata);

CurForLength  = ADIGATORFORDATA(CSFloc).FOR(1).LENGTHS(ForIter);
CurDataLength = NUMCSF/CurForLength;
if NUMCSF && CurForLength && ceil(CurDataLength) ~= CurDataLength
  % For some reason the data is short. This is probably due to an IF
  % statement that didnt fire on the last iteration of the loop.
  CurDataLength = ceil(CurDataLength);
  NUMCSF        = CurDataLength*CurForLength;
  CSFdata(NUMCSF,1) = 0;
end
if ForIter == 1
  % First Iteration of this Loop.
  ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr) = CSFdata;
elseif ~NUMCSF
  % This Special Function never got called (either in an IF statement or a
  % nested FOR loop)
  if ~isempty(OldData)
    % This function got called on a previous iteration, use 0 as a
    % placeholder
    ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)(end,ForIter) = 0;
  end
elseif isempty(OldData)
  % First iteration which this function fired - just store it in proper
  % iteration column
  ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)(1:NUMCSF,ForIter)...
    = CSFdata;
else
  % Check For child FORLENGTH and/or indice sizes changing
  if CurForLength == 0
    % This FOR loop never fired
    ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)(:,ForIter) = 0;
    return
  end

  PrevForLength = max(ADIGATORFORDATA(CSFloc).FOR(1).LENGTHS(1:ForIter-1));
  if PrevForLength == 0
    ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)...
      (1:NUMCSF,ForIter) = CSFdata;
    return
  end
  PrevDataLength = size(ADIGATORFORDATA(ForCount).(FunStr)...
    (SFcount).(RefStr),1)/PrevForLength;
  if CurForLength ~= PrevForLength || CurDataLength ~= PrevDataLength
    NewInds = ReshapeData(CurForLength,PrevForLength,CurDataLength,...
      PrevDataLength,CSFdata,ForIter,OldData);
    ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)...
      (1:size(NewInds,1),1:size(NewInds,2)) = NewInds;
  else
    ADIGATORFORDATA(ForCount).(FunStr)(SFcount).(RefStr)...
      (1:NUMCSF,ForIter) = CSFdata;
  end
end
% Empty out the Data spot.
ADIGATORFORDATA(CSFloc).(FunStr)(CSFcount).(RefStr) = [];
return
end

function NewData = ReshapeData(CurForLength,PrevForLength,...
  CurDataLength,PrevDataLength,CSFdata,ForIter,OldData)
% Reshapes Data - used for instances when either FOR lengths are changing
% or indices are changing size.

if CurForLength > PrevForLength
  % Current Child FORLENGTH is greater than Previous ones
  if CurDataLength > PrevDataLength
    % Current Child data length is greater than previous
    NewData = zeros(CurDataLength*CurForLength,ForIter);
    for Icount = 1:ForIter-1
      TempData = OldData(:,Icount);
      TempData = reshape(TempData,PrevDataLength,PrevForLength);
      TempData(CurDataLength,CurForLength) = 0;
      NewData(:,Icount) = TempData(:);
    end
    NewData(:,ForIter) = CSFdata;
  else
    % Current Child data length is <= previous
    NewData = zeros(PrevDataLength*CurForLength,ForIter);
    for Icount = 1:ForIter-1
      TempData = OldData(:,Icount);
      TempData = reshape(TempData,PrevDataLength,PrevForLength);
      TempData(PrevDataLength,CurForLength) = 0;
      NewData(:,Icount) = TempData(:);
    end
    if CurDataLength < PrevDataLength
      % Current Child data length is < previous
      TempData = reshape(CSFdata,CurDataLength,CurForLength);
      TempData(PrevDataLength,CurForLength) = 0;
      NewData(:,ForIter) = TempData(:);
    else
      NewData(:,ForIter) = CSFdata(:);
    end
  end
elseif CurForLength < PrevForLength
  % Current Child FORLENGTH is less than previous
  if CurDataLength > PrevDataLength
    NewData = zeros(PrevForLength*CurDataLength,ForIter);
    for Icount = 1:ForIter-1
      TempData = OldData(:,Icount);
      TempData = reshape(TempData,PrevDataLength,PrevForLength);
      TempData(CurDataLength,PrevForLength) = 0;
      NewData(:,Icount) = TempData(:);
    end
    NewData(1:length(CSFdata),ForIter) = CSFdata;
  else
    % Old Data stays the same
    CSFdata = reshape(CSFdata,CurDataLength,CurForLength);
    CSFdata(PrevDataLength,PrevForLength) = 0;
    NewData = [OldData CSFdata(:)];
  end
else
  % FORLENGTHS are the same
  if CurDataLength > PrevDataLength
    NewData = zeros(CurDataLength*CurForLength,ForIter);
    for Icount = 1:ForIter-1
      TempData = OldData(:,Icount);
      TempData = reshape(TempData,PrevDataLength,PrevForLength);
      TempData(CurDataLength,PrevForLength) = 0;
      NewData(:,Icount) = TempData(:);
    end
    NewData(:,ForIter) = CSFdata(:);
  else
    CSFdata = reshape(CSFdata,CurDataLength,CurForLength);
    CSFdata(PrevDataLength,CurForLength) = 0;
    NewData = [OldData CSFdata(:)];
  end
end
% RefLength CRinds Inds
return
end

function [outEvalStr, outEvalVar] = DoRemapping(VarCounts,ForCount)
global ADIGATOR ADIGATORVARIABLESTORAGE
SaveCounts = VarCounts(logical(ADIGATOR.VARINFO.SAVE.FOR(VarCounts,2)));

nOutEval   = length(SaveCounts);
outEvalStr = cell(nOutEval,1);
outEvalVar = cell(nOutEval,1);
for Vcount = 1:nOutEval
  Scount  = SaveCounts(Vcount);
  NameLoc = ADIGATOR.VARINFO.NAMELOCS(Scount,1);
  VarStr  = ADIGATOR.VARINFO.NAMES{NameLoc};
  
  OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(Scount,1);
  OverVar = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
  SaveLoc = ADIGATOR.VARINFO.SAVE.FOR(Scount,2);
  SaveVar = ADIGATORVARIABLESTORAGE.SAVE{SaveLoc};
  
  ForSaveLoc = ADIGATOR.VARINFO.SAVE.FOR(Scount,1);
  % Check for Loop inside of Conditional Statement
  if ~isempty(ADIGATOR.VARINFO.OVERMAP.IF)
    IfOverLoc = ADIGATOR.VARINFO.OVERMAP.IF(Scount,1);
    IfLoopLoc = ADIGATOR.VARINFO.OVERMAP.IF(Scount,2);
    IfSaveLoc = ADIGATOR.VARINFO.SAVE.IF(Scount,1);
    if IfLoopLoc == ForCount
      % Loop Inside of Conditonal Statement
       CondOvermap = ADIGATORVARIABLESTORAGE.OVERMAP{IfOverLoc};
       if ~ADIGATOR.VARINFO.SAVE.IF(Scount,2)
        % Can do the remapping here
        CondOvermap = cadaPrintReMap(OverVar,CondOvermap,Scount);
        % Return the Conditional Overmap
        outEvalVar{Vcount} = CondOvermap;
       else
         % Need to remap to the end of the loop, and store the end of the
         % loop
         if isa(OverVar,'cada')
           outEvalVar{Vcount} = cadaPrintReMap(OverVar,SaveVar,Scount);
         end
         ADIGATORVARIABLESTORAGE.SAVE{ADIGATOR.VARINFO.SAVE(Scount,2)} = ...
           outEvalVar{Vcount};
       end
      
       if IfSaveLoc
         ADIGATORVARIABLESTORAGE.SAVE{IfSaveLoc} = CondOvermap;
       end
    else
      % Conditional Statement Inside of Loop or Not a Conditional Statement
      % at all
      if isa(OverVar,'cada')
        outEvalVar{Vcount} = cadaPrintReMap(OverVar,SaveVar,Scount);
      end
      if IfSaveLoc
        ADIGATORVARIABLESTORAGE.SAVE{IfSaveLoc} = outEvalVar{Vcount};
      end
    end
  else
    outEvalVar{Vcount} = cadaPrintReMap(OverVar,SaveVar,Scount);
  end
  if ForSaveLoc
    ADIGATORVARIABLESTORAGE.SAVE{ForSaveLoc} = outEvalVar{Vcount};
  end
  outEvalStr{Vcount} = sprintf([VarStr,' = adigatorForEvalVar{%1.0d};'],Vcount);
end

end

function [outEvalStr, outEvalVar] = getContBreakOvermaps(OverLocs,CBflag,LoopCounts)
global ADIGATOR ADIGATORVARIABLESTORAGE

outEvalVar = ADIGATORVARIABLESTORAGE.OVERMAP(OverLocs);
nOutEval   = length(OverLocs);
outEvalStr = cell(nOutEval,1);
if ~CBflag
  PrevCounts = 1:LoopCounts(1)-1;
  PostCounts = LoopCounts(end)+1:size(ADIGATOR.VARINFO.OVERMAP.IF,1);
end
for Vcount = 1:nOutEval
  OverLoc = OverLocs(Vcount);
  if CBflag
    VarCount = LoopCounts(find(ADIGATOR.VARINFO.OVERMAP.CONT(LoopCounts) == OverLoc,1,'last'));
  else
    VarCount = LoopCounts(find(ADIGATOR.VARINFO.OVERMAP.BREAK(LoopCounts) == OverLoc,1,'last'));
  end
  if ~CBflag
    % BREAK STATEMENTS
    IfOverLoc = ADIGATOR.VARINFO.OVERMAP.IF(VarCount,1);
    if IfOverLoc && (...
        any(ADIGATOR.VARINFO.OVERMAP.IF([PrevCounts, PostCounts],1)==IfOverLoc) ||...
        any(ADIGATOR.VARINFO.OVERMAP.IF(PrevCounts,2)==IfOverLoc) )
      % This variable belongs to an overmap for a conditional statement
      % outside of the loop - assign that overmap.
      IfOverMap = ADIGATORVARIABLESTORAGE.OVERMAP{IfOverLoc};
      BreakOverMap = ADIGATORVARIABLESTORAGE.OVERAMP{OverLoc};
      Union = cadaUnionVars(IfOverMap,BreakOverMap);
      ADIGATORVARIABLESTORAGE.OVERMAP{IfOverLoc} = Union;
      ReturnLoc = ADIGATOR.VARINFO.RETURN(VarCount);
      if ReturnLoc < 0
        % Return a previously saved variable
        outEvalVar{Vcount} = ADIGATORVARIABLESTORAGE.SAVE{-ReturnLoc};
      elseif ReturnLoc
        % Return the CONDITIONAL OVERMAP
        outEvalVar{Vcount} = Union;
      end
    else
      Union = [];
    end
    SaveLoc = ADIGATOR.VARINFO.SAVE.FOR(VarCount,2);
    if SaveLoc
      % This variable is set to be saved so that we can re-map to it during
      % the printing evaluations - we want to save this to be the BREAK
      % OVERMAP
      if ~isempty(Union)
        % This belongs to a conditional overmap from the case above, use
        % that variable
        ADIGATORVARIABLESTORAGE.SAVE{SaveLoc} = Union;
      else
        ADIGATORVARIABLESTORAGE.SAVE{SaveLoc} = outEvalVar{Vcount};
      end
    end
  end
  
  NameLoc = ADIGATOR.VARINFO.NAMELOCS(VarCount,1);
  VarStr  = ADIGATOR.VARINFO.NAMES{NameLoc};
  outEvalStr{Vcount} = sprintf([VarStr,' = adigatorForEvalVar{%1.0d};'],Vcount);
end

end
