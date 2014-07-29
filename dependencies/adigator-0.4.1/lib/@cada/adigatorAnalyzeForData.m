function adigatorAnalyzeForData(FORCOUNT,dummyVar) %#ok<INUSD>
%function adigatorAnalyzeForData(FORCOUNT,StartVars)
% This routine is called from either the main adigator file (in this case it
% is being called for a user sub-function) or the
% adigatorFlowControlAnalyzer0 (called for an outer FOR loop in the main
% user function). This routine sorts through the data obtained in the
% Pre-Printing Evaluation for each organizational function to determine
% which of these will require special function and/or derivative
% reference/assignment indexing. For the most part, for any organizational
% function, if either the sizes of the inputs/outputs or the
% reference/assignment indices (if applicable) change size, then some
% additional information must be printed.
% ----------------------- Input Information ----------------------------- %
% FORCOUNT:   the element of the ADIGATORFORDATA structure array we are
%             looking at
% dummyVar:   just a dummy variable so that we can store this routine
%             within the overloaded folder.
% -------------------- Sub-Routine Listing ------------------------------ %
% GetDataDependencies: given a set of indices and/or sizes collected over
%    multiple embedded FOR loops, this routine recursively calls itself to
%    find which loops the data is dependent upon
% RemoveUnnneededIndices: given the dependency map from
%    GetDataDependencies, this routine removes the ones which are repeated
%    (i.e. independent)
% GetSubsrefInds:   obtains the derivative index mapping for subsref
% GetSubsasgnInds:  obtains the derivative index mapping for subsasgn
% GetHorzcatInds:   obtains the derivative index mapping for horzcat
% GetVertcatInds:   obtains the derivative index mapping for vertcat
% GetTransposeInds: obtains the derivative index mapping for transpose
% GetRepmatInds:    obtains the derivative index mapping for repmat
% GetReshapeInds:   obtains the derivative index mapping for reshape
% AssignReferenceNames: used to assign reference names to the
%   ADIGATORFORDATA structure
% AssignReferenceNamesCat: version of AssignReferenceNames only used with
% horzcat and vertcat (needed since they have an extra unknown dimension)
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0


global ADIGATOR ADIGATORDATA ADIGATORFORDATA ADIGATORVARIABLESTORAGE
DataNameCount = 0;
NUMvod        = ADIGATOR.NVAROFDIFF;
NDstr         = sprintf('%1.0d',ADIGATOR.DERNUMBER);

% ----------------------------------------------------------------------- %
%                EMBEDDED FOR LOOP ITERATION SIZES                        %
% ----------------------------------------------------------------------- %
ADIGATORFORDATA(FORCOUNT).MAXLENGTH = ADIGATORFORDATA(FORCOUNT).FOR(1).LENGTHS;
OuterLoopMaxLength = ADIGATORFORDATA(FORCOUNT).MAXLENGTH;
for Fcount = 2:length(ADIGATORFORDATA(FORCOUNT).FOR)
  % -------------------- Check for Dependencies ------------------------- %
  ForLengths  = ADIGATORFORDATA(FORCOUNT).FOR(Fcount).LENGTHS;
  FORLOCS     = ADIGATORFORDATA(FORCOUNT).FOR(Fcount).LOCS;
  if size(ForLengths,2) < OuterLoopMaxLength
    ForLengths(end,OuterLoopMaxLenght) = 0;
  end
  % Set MAXLENGTH field of ADIGATORFORDATA
  ADIGATORFORDATA(FORLOCS(1,end)).MAXLENGTH = max(ForLengths(:));
  NUMloops    = size(FORLOCS,2)-1;
  FORLENGTHS  = zeros(1,NUMloops);
  for F2count = 1:NUMloops
    FORLENGTHS(F2count) = ADIGATORFORDATA(FORLOCS(1,F2count)).MAXLENGTH;
  end
  FORDEP = GetDataDependencies(ForLengths,FORLENGTHS);
  % .LENGTHS field will now be a 1x2 cell array - set it as such
  for F2count = 1:size(FORLOCS,2)
    ADIGATORFORDATA(FORLOCS(1,F2count)).FOR(FORLOCS(2,F2count)).LENGTHS...
      = cell(1,2);
  end
  if nnz(FORDEP) > 0
    % --------- Is Dependent - Assign Printed Values/Names -------------- %
    DepLengths = RemoveUnneededData(ForLengths,FORLENGTHS,FORDEP);

    % If an element of ADIGATORFORDATA has a FOR.LENGTHS field with strings
    % in FOR(i>1).LENGTHS{1} and FOR(i>1).LENGTHS{2}, then it prints:
    %       FOR(i).LENGTHS{1} = FOR(i).LENGTHS{2}; this will be a reshape
    %       or a reference.
    % If an element of ADIGATORFORDATA has a FOR.LENGTHS field with a string
    % in FOR(1).LENGTHS{1}, then it uses this to reference off of its
    % looping variable to tell the code that the iterations change size.
    DepForLengths = FORLENGTHS(FORDEP);
    DepForLocs    = FORLOCS(1,FORDEP);
    if FORDEP(1); ADIGATOR.FUNDEP = 1; end
    
    if ~isempty(DepForLocs)
      Inds = reshape(DepLengths,numel(DepLengths)/DepForLengths(1),DepForLengths(1));
      CountName = ADIGATORFORDATA(DepForLocs(1)).COUNTNAME;
      DepForLocs(1) = []; DepForLengths(1) = [];
    end
    
    % ------------------------ Assign Data Name ------------------------- %
    if nnz(Inds) < numel(Inds); Inds = sparse(Inds); end
    INDEXCOUNT = ADIGATORDATA.INDEXCOUNT + 1;
    ADIGATORDATA.INDEXCOUNT = INDEXCOUNT;
    ADIGATORDATA.INDICES.(sprintf('Index%1.0d',INDEXCOUNT)) = Inds;
    IndName = sprintf('Gator%1.0dIndices.Index%1.0d',...
      ADIGATOR.DERNUMBER,INDEXCOUNT);

    % Middle FOR loops - reshaping/referencing stuff
    for F2count = 1:length(FORLENGTHS)
      if FORDEP(F2count)
        if ~isempty(DepForLengths)
          Inds = reshape(Inds(:,1),numel(Inds(:,1))/DepForLengths(1),DepForLengths(1));
          DepForLengths(1) = [];
          ADIGATORFORDATA(FORLOCS(1,F2count+1)).FOR(FORLOCS(2,F2count+1)).LENGTHS{2} = ...
            sprintf(['reshape(',IndName,'(:,',CountName,'),%1.0f,%1.0f)']...
            ,size(Inds,1),size(Inds,2));
          CountName = ADIGATORFORDATA(DepForLocs(1)).COUNTNAME;
          DepForLocs(1) = [];
        else
          ADIGATORFORDATA(FORLOCS(1,F2count+1)).FOR(FORLOCS(2,F2count+1)).LENGTHS{2} = ...
            [IndName,'(',CountName,')'];
        end
        DataNameCount = DataNameCount+1;
        IndName = sprintf(['cada',NDstr,'forindex%1.0f'],DataNameCount);
        ADIGATORFORDATA(FORLOCS(1,F2count+1)).FOR(FORLOCS(2,F2count+1)).LENGTHS{1} = IndName;
      end
    end
    % Innermost FOR loop - what adigatorFlowControlAnalyzer2 will see when
    % printing the FOR loop.
    ADIGATORFORDATA(FORLOCS(1,end)).FOR(FORLOCS(2,end)).LENGTHS{1}...
      = IndName;
  end
end

% ----------------------------------------------------------------------- %
%                              SUBSREF                                    %
% ----------------------------------------------------------------------- %
for Rcount   = 1:length(ADIGATORFORDATA(FORCOUNT).SUBSREF)
  % ---- Check for Input Size and/or Reference Index Dependencies ------- %
  % We use the syntax y = x(inds);
  % Indices have been collected throughout - sizes are stacked upon one
  % another from within the innermost loop.
  REFINDS    = ADIGATORFORDATA(FORCOUNT).SUBSREF(Rcount).INDICES;
  REFLOCS    = ADIGATORFORDATA(FORCOUNT).SUBSREF(Rcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,REFLOCS);
  if isempty(REFINDS)
    % This call to subsref never fires for some reason
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(REFLOCS(1,Fcount)).SUBSREF(REFLOCS(2,Fcount))...
        .INDICES = [];
      ADIGATORFORDATA(REFLOCS(1,Fcount)).SUBSREF(REFLOCS(2,Fcount))...
        .SIZES   = [];
    end
    continue
  elseif size(REFINDS,2) < OuterLoopMaxLength
    REFINDS(end,OuterLoopMaxLength) = 0;
  end
  
  ySizes = ADIGATORFORDATA(REFLOCS(1,end)).SUBSREF(REFLOCS(2,end)).SIZES...
    (1:2,:);
  xSizes = ADIGATORFORDATA(REFLOCS(1,end)).SUBSREF(REFLOCS(2,end)).SIZES...
    (3:4,:);

  % Pre-Allocate SUBSREF .INDICES and .SIZES - We use .INDICES field to
  % store information on derivative indexing and .SIZES field to store
  % information on variables changing size.
  for Fcount = 1:length(FORLENGTHS)
    ADIGATORFORDATA(REFLOCS(1,Fcount)).SUBSREF(REFLOCS(2,Fcount))...
      .INDICES = cell(NUMvod,3);
    ADIGATORFORDATA(REFLOCS(1,Fcount)).SUBSREF(REFLOCS(2,Fcount))...
      .SIZES   = cell(3,3);
  end
  
  % get overmapped x and y (y = x(inds))
  xOver = ADIGATORFORDATA(REFLOCS(1,end)).SUBSREF(REFLOCS(2,end)).VARS{1};
  if ~isa(xOver,'cada')
    if isempty(xOver)
      % xLoc will be empty whenever we are in a 2nd derivative and we are
      % referencing off of a set of FOR loop indices - dont worry about it,
      % it is just numeric references and none of the references ever change
      % size.
      continue
    else
      xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
    end
  end
  yOver = ADIGATORFORDATA(REFLOCS(1,end)).SUBSREF(REFLOCS(2,end)).VARS{2};
  
  % Maximum number of references per SUBSREF call
  NUMinds = numel(REFINDS)/prod(FORLENGTHS);
  
  %                   Function Sizing Checks                              %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  % Need to check for two things:
  %     1. X rows changing size, where X is a matrix. - Need to
  %     reorder the reference indices, so proper derivatives are
  %     found. Also, if there is a Linear Index, need to correct
  %     this in the FUNCTION reference.
  %     2. Y changing size on either dimension - need to print
  %     stuff out so that the proper reference can be done on the
  %     FUNCTION, does not affect the derivative.
  
  % Check to see if xCols > 1 AND xRows changes size
  XSIZES = [];
  XSIZEDEP = [];
  if xOver.func.size(2) > 1 && ~any(isinf(xOver.func.size)) ...
      && any(xSizes(1,:) ~= xOver.func.size(1))
    % Indexing is no longer proper for cases where xRowSize <
    % max(xRowSize) - need to change it
    NewInds = reshape(REFINDS,NUMinds,prod(FORLENGTHS));
    RScount = 0; MaxRef = zeros(xOver.func.size(1),xOver.func.size(2));
    MaxRef(:) = 1:numel(MaxRef);
    for Icount = 1:size(NewInds,2)
      if sum(NewInds(:,Icount))
        RScount = RScount+1;
        xRef = MaxRef(1:xSizes(1,RScount),1:xSizes(2,RScount));
        OldInd = NewInds(:,Icount);
        NewInds(logical(OldInd),Icount) = xRef(nonzeros(OldInd));
      end
    end
    REFINDS = reshape(NewInds,size(REFINDS,1),size(REFINDS,2));
    if ~ADIGATORFORDATA(REFLOCS(1,end)).SUBSREF(REFLOCS(2,end)).FLAGS(2)
      % Linear indexing is used - will need to print out more stuff
      % to fix this in the output file so that the FUNCTION is
      % referenced properly. If a linear index is not used, don't
      % need to.
      xSize = xSizes(1,:);
      % First dimension changes - put it into same form as indice
      % array.
      xSizeTemp = zeros(1,prod(FORLENGTHS));
      xSizeTemp(logical(sum(reshape(REFINDS,NUMinds,prod(FORLENGTHS)),1)))...
        = xSize;
      xSize = reshape(xSizeTemp,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % Get SIZE dependencies
      XSIZEDEP = GetDataDependencies(xSize,FORLENGTHS);
      % Remove any unneeded sizes
      XSIZES = RemoveUnneededData(xSize,FORLENGTHS,XSIZEDEP);
    end
  end
  

  % Check for changing y-size - used to make the proper FUNCTION reference
  % in the overloaded function.
  YSIZES = cell(2,1);
  YSIZEDEP = cell(2,1);
  for RScount = 1:2
    ySize = ySizes(RScount,:);
    if sum(ySize~=ySize(1))
      % This dimension changes - put it into same form as indice
      % array.
      ySizeTemp = zeros(1,prod(FORLENGTHS));
      ySizeTemp(logical(sum(reshape(REFINDS,NUMinds,prod(FORLENGTHS)),1)))...
        = nonzeros(ySize);
      ySize = reshape(ySizeTemp,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % Get SIZE dependencies
      YSIZEDEP{RScount} = GetDataDependencies(ySize,FORLENGTHS);
      % Remove any unneeded sizes
      YSIZES{RScount} = RemoveUnneededData(ySize,FORLENGTHS,YSIZEDEP{RScount});
    end
  end
  
  %                    Get the Derivative Indices                         %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  
  % Get SUBSREF indice Dependencies
  REFDEP = GetDataDependencies(REFINDS,FORLENGTHS);
  % Remove any unneeded indices
  DEPINDS = RemoveUnneededData(REFINDS,FORLENGTHS,REFDEP);
  % Check for b changing size and bCols > 1
  if yOver.func.size(2) > 1 && any(ySizes(1,:)~= yOver.func.size(1)) && ...
      ~any(isinf(yOver.func.size))
    % y changes row size and is a matrix, need to put ySizes into
    % same form as DEPINDS so that we can properly get the
    % derivative indices.
    ySizeTemp = zeros(2,prod(FORLENGTHS));
    ySizeTemp(:,logical(sum(reshape(REFINDS,NUMinds,prod(FORLENGTHS)),1)))...
      = ySizes;
    ySizeTemp = reshape(ySizeTemp,2*prod(FORLENGTHS(2:end)),FORLENGTHS(1));
    ySizes = RemoveUnneededData(ySizeTemp,FORLENGTHS,REFDEP);
    if size(ySizes,2) > 1
      ySizes = reshape(ySizes,2,prod(FORLENGTHS(REFDEP)));
    end
  else
    ySizes = [];
  end
  
  % Get the Derivative Indices.
  if size(DEPINDS,2) > 1
    DEPINDS = reshape(DEPINDS,NUMinds,prod(FORLENGTHS(REFDEP)));
  end
  
  DEPINDS(DEPINDS < 0) = 0;
  YSIZES{1}(YSIZES{1} < 0) = 0;
  YSIZES{2}(YSIZES{2} < 0) = 0;
  ySizes(ySizes < 0) = 0;

  DerInds = GetSubsrefInds(xOver,yOver,DEPINDS,ySizes); 
  %                Assigning Names for Printing Run                       %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  
  % Assign Derivative Indice names that will get printed
  for Vcount = 1:NUMvod
    if ~isempty(DerInds{Vcount})
      % These get assigned to .INDICES{Vcount,1:3}
      DataNameCount = AssignReferenceNames(DerInds{Vcount},REFDEP,REFLOCS,...
        FORLENGTHS,'SUBSREF','INDICES',Vcount,0,DataNameCount);
    end
  end
  % Assign Function Size names that will get printed
  if ~isempty(XSIZES)
    % X First Dimension Changing Size - print it out. 
    % These get assigned to .SIZES{1,1:3}
    DataNameCount = AssignReferenceNames(XSIZES,XSIZEDEP,REFLOCS,...
      FORLENGTHS,'SUBSREF','SIZES',1,0,DataNameCount);
  end
  if ~isempty(YSIZES{1})
    % Y First Dimension Changing Size - print it out.
    % These get assigned to .SIZES{2,1:3}
    DataNameCount = AssignReferenceNames(YSIZES{1},YSIZEDEP{1},REFLOCS,...
      FORLENGTHS,'SUBSREF','SIZES',2,0,DataNameCount);
  end
  if ~isempty(YSIZES{2})
    % Y Second Dimension Changing Size - print it out.
    % These get assigned to .SIZES{3,1:3}
    DataNameCount = AssignReferenceNames(YSIZES{2},YSIZEDEP{2},REFLOCS,...
      FORLENGTHS,'SUBSREF','SIZES',3,0,DataNameCount);
  end
end

% ----------------------------------------------------------------------- %
%                              SUBSASGN                                   %
% ----------------------------------------------------------------------- %
for Scount   = 1:length(ADIGATORFORDATA(FORCOUNT).SUBSASGN)
  % Syntax: X(inds) = B;
  % Indices have been collected throughout - sizes are stacked upon one
  % another from within the innermost loop.
  ASGNINDS   = ADIGATORFORDATA(FORCOUNT).SUBSASGN(Scount).INDICES;
  ASGNLOCS   = ADIGATORFORDATA(FORCOUNT).SUBSASGN(Scount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,ASGNLOCS);
  if isempty(ASGNINDS)
    % This call to subsasgn never fires
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SUBSASGN(ASGNLOCS(2,Fcount))...
        .INDICES = [];
      ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SUBSASGN(ASGNLOCS(2,Fcount))...
        .SIZES   = [];
    end
    continue
  elseif size(ASGNINDS,2) < OuterLoopMaxLength
    ASGNINDS(end,OuterLoopMaxLength) = 0;
  end
  
  % Sizes of X for each loop iteration
  xSizes   = ...
    ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).SIZES(1:2,:);
  % Sizes of B for each loop iteration
  bSizes   = ...
    ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).SIZES(3:4,:);
  % get overmapped x and b
  xOverLoc = ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).VARS{1};
  xOver    = ADIGATORVARIABLESTORAGE.OVERMAP{xOverLoc};

  bOver    = ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).VARS{2};
  if ~isa(bOver,'cada')
    bOver  = ADIGATORVARIABLESTORAGE.OVERMAP{bOver};
    ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).VARS{2} = bOver;
  end
  
  % Maximum number of references per SUBSASGN call
  NUMinds = numel(ASGNINDS)/prod(FORLENGTHS);
  
  %                    Function Sizing Checks                             %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  % Need to Check for two things:
  %     1. X changing Row size, when X is a matrix. In this case
  %     need to re-order the assignment indices so that they match
  %     X's OverMapped Linear Index. Also, if a linear indexing is
  %     used in the assignment, need to change this in the FUNCTION
  %     reference.
  %     2. B changing Row or Column size, need to print some stuff
  %     out to ensure proper FUNCTION assignment in the output
  %     file. This does not affect the derivatives.
  
  XSIZES = [];
  XSIZEDEP = [];
  % Check for xRow varying and xCol > 1
  if xOver.func.size(2) > 1 && any(xSizes(1,:)~=xOver.func.size(1)) && ...
      ~any(xOver.func.size)
    % the indices in ASGNINDS are no longer valid, need to convert
    % them to match the overmapped X.
    NewInds = reshape(ASGNINDS,NUMinds,prod(FORLENGTHS));
    SScount = 0; MaxRef = zeros(xOver.func.size(1),xOver.func.size(2));
    MaxRef(:) = 1:numel(MaxRef);
    for Icount = 1:size(NewInds,2)
      if sum(NewInds(:,Icount))
        SScount = SScount+1;
        xRef = MaxRef(1:xSizes(1,SScount),1:xSizes(2,RScount));
        OldInd = NewInds(:,Icount);
        NewInds(logical(OldInd),Icount) = xRef(nonzeros(OldInd));
      end
    end
    ASGNINDS = reshape(NewInds,size(ASGNINDS,1),size(ASGNINDS,2));
    % ASGNINDS are now changed to match OverMap.
    if ~ADIGATORFORDATA(ASGNLOCS(1,end)).SUBSASGN(ASGNLOCS(2,end)).FLAGS(2)
      % Linear indexing is used - will need to print out more stuff
      % to fix this in the output file so that the FUNCTION is
      % assigned properly. If a linear index is not used, don't
      % need to.
      xSize = xSizes(1,:);
      % First dimension changes - put it into same form as indice
      % array.
      xSizeTemp = zeros(1,prod(FORLENGTHS));
      xSizeTemp(logical(sum(reshape(ASGNINDS,NUMinds,prod(FORLENGTHS)),1)))...
        = xSize;
      xSize = reshape(xSizeTemp,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % Get SIZE dependencies
      XSIZEDEP = GetDataDependencies(xSize,FORLENGTHS);
      % Remove any unneeded sizes
      XSIZES = RemoveUnneededData(xSize,FORLENGTHS,XSIZEDEP);
    end
  end

  BSIZES = cell(2,1);
  BSIZEDEP = cell(2,1);
  % Check for changing b-size
  for SScount = 1:2
    bSize = bSizes(SScount,:);
    if any(bSize ~= bSize(1))
      % This dimension changes - put it into same form as indice
      % array.
      % Zeros in bSize mean there was an empty assignment, make these into
      % Inf so that GetDataDependencies and RemoveUnneededIndices doesnt
      % think that they are cases which subsasgn just didnt fire.
      bSizeTemp = zeros(1,prod(FORLENGTHS));
      bSizeTemp(logical(sum(reshape(ASGNINDS,NUMinds,prod(FORLENGTHS)),1)))...
        = nonzeros(bSize);
      bSize = reshape(bSizeTemp,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % Get SIZE dependencies
      BSIZEDEP{SScount} = GetDataDependencies(bSize,FORLENGTHS);
      % Remove any unneeded sizes
      BSIZES{SScount} = RemoveUnneededData(bSize,FORLENGTHS,...
        BSIZEDEP{SScount});
      BSIZES{SScount}(BSIZES{SScount}==inf) = 0;
    end
  end
  %                    Get the Derivative Indices                         %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  % This gets all of the derivative indices which need to get
  % printed. - GetSubsasgnInds does all the hard work.
  % Get SUBSASGN indice Dependencies
  [ASGNDEP] = GetDataDependencies(ASGNINDS,FORLENGTHS);
  % Remove any unneeded indices
  DEPINDS = RemoveUnneededData(ASGNINDS,FORLENGTHS,ASGNDEP);
  DEPINDS(DEPINDS==inf) = 0;

  % Check for b changing size and bCols > 1
  if bOver.func.size(2) > 1 && any(bSizes(1,:) ~= bOver.func.size(1)) && ...
      ~any(isinf(bOver.func.size))
    % b changes row size and is a matrix, need to put bSizes into
    % same form as DEPINDS so that we can properly get the
    % derivative indices.
    bSizeTemp = zeros(2,prod(FORLENGTHS));
    bSizeTemp(:,logical(sum(reshape(ASGNINDS,NUMinds,prod(FORLENGTHS)),1))) = bSizes;
    bSizeTemp = reshape(bSizeTemp,2*prod(FORLENGTHS(2:end)),FORLENGTHS(1));
    bSizes = RemoveUnneededData(bSizeTemp,FORLENGTHS,ASGNDEP);
    if size(bSizes,2) > 1
      bSizes = reshape(bSizes,2,prod(FORLENGTHS(ASGNDEP)));
    end
  else
    bSizes = [];
  end
  % Get the Derivative indices
  if size(DEPINDS,2) > 1
    DEPINDS = reshape(DEPINDS,NUMinds,prod(FORLENGTHS(ASGNDEP)));
  end
  DerInds = GetSubsasgnInds(xOver,bOver,DEPINDS,bSizes);
  
  
  %                   Assign Names for Printing Run                       %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  
  % Pre-Allocate SUBSASGN .INDICES and .SIZES fields - We use the .INDICES field to
  % store derivative indexing information and .SIZES field to store
  % information on inputs/outputs changing size
  for Fcount = 1:length(FORLENGTHS)
    ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SUBSASGN(ASGNLOCS(2,Fcount))...
      .INDICES = cell(NUMvod,6);
    ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SUBSASGN(ASGNLOCS(2,Fcount))...
      .SIZES   = cell(3,3);
  end
  
  % Assign Derivative Indice names that will get printed
  for Vcount = 1:NUMvod
    if ~isempty(DerInds{Vcount,1})
      % Assignment Indices - These get stored in .INDICES{Vcount,1:3}
      DataNameCount = AssignReferenceNames(DerInds{Vcount,1},ASGNDEP,...
        ASGNLOCS,FORLENGTHS,'SUBSASGN','INDICES',Vcount,0,DataNameCount);
    end    
    if ~isempty(DerInds{Vcount,2}) && sum(DerInds{Vcount,2}(:))
      % Derivative Removal Indices - These get stored in
      % .INDICES{Vcount,4:6}
      DataNameCount = AssignReferenceNames(DerInds{Vcount,2},ASGNDEP,...
        ASGNLOCS,FORLENGTHS,'SUBSASGN','INDICES',Vcount,3,DataNameCount);
    end
  end
  
  % Assign Function Size names that will get printed
  if ~isempty(XSIZES)
    % X First Dimension Changing Size - print it out. 
    % These get stored in .SIZES{1,1:3}
    DataNameCount = AssignReferenceNames(XSIZES,XSIZEDEP,ASGNLOCS,...
      FORLENGTHS,'SUBSASGN','SIZES',1,0,DataNameCount);
  end
  if ~isempty(BSIZES{1})
    % B First Dimension Changing Size - print it out.
    % These get stored in .SIZES{2,1:3}
    DataNameCount = AssignReferenceNames(BSIZES{1},BSIZEDEP{1},ASGNLOCS,...
      FORLENGTHS,'SUBSASGN','SIZES',2,0,DataNameCount);
  end
  if ~isempty(BSIZES{2})
    % Y Second Dimension Changing Size - print it out.
    % These get stored in .SIZES{3,1:3}
    DataNameCount = AssignReferenceNames(BSIZES{2},BSIZEDEP{2},ASGNLOCS,...
      FORLENGTHS,'SUBSASGN','SIZES',3,0,DataNameCount);
  end
end
% ----------------------------------------------------------------------- %
%                              SPARSE                                     %
% ----------------------------------------------------------------------- %
for Scount   = 1:length(ADIGATORFORDATA(FORCOUNT).SPARSE)
  % y = sparse(rows,cols,x,M,N) - This is just a SubsAsgn into a
  % zero matrix.
  ASGNINDS   = ADIGATORFORDATA(FORCOUNT).SPARSE(Scount).INDICES;
  ASGNLOCS   = ADIGATORFORDATA(FORCOUNT).SPARSE(Scount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,ASGNLOCS);
  
  if isempty(ASGNINDS)
    % This call to SPARSE never fires
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SPARSE(ASGNLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,3);
      ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SPARSE(ASGNLOCS(2,Fcount))...
        .SIZES   = cell(2,3);
    end
    continue
  elseif size(ASGNINDS,2) < OuterLoopMaxLength
    ASGNINDS(end,OuterLoopMaxLength) = 0;
  end
  % Sizes of Y for each loop iteration
  ySizes = ADIGATORFORDATA(ASGNLOCS(1,end)).SPARSE(ASGNLOCS(2,end)).SIZES(1:2,:);
  % Sizes of X for each loop iteration
  xSizes = ADIGATORFORDATA(ASGNLOCS(1,end)).SPARSE(ASGNLOCS(2,end)).SIZES(3:4,:);
  % get overmapped Y and X
  yOver = ADIGATORFORDATA(ASGNLOCS(1,end)).SPARSE(ASGNLOCS(2,end)).VARS{1};
  xOver = ADIGATORFORDATA(ASGNLOCS(1,end)).SPARSE(ASGNLOCS(2,end)).VARS{2};
  if ~isa(xOver,'cada')
    xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
    ADIGATORFORDATA(ASGNLOCS(1,end)).SPARSE(ASGNLOCS(2,end)).VARS{2} = xOver;
  end
  
  % Maximum number of references per SPARSE call
  NUMinds = numel(ASGNINDS)/prod(FORLENGTHS);
  
  %                    Function Sizing Checks                             %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  % Need to Check for two things:
  %     1. Y changing Row size, when Y is a matrix. In this case
  %     need to re-order the assignment indices so that they match
  %     Y's OverMapped Linear Index.
  %     2. X changing Row or Column size, need to print some stuff
  %     out to ensure proper FUNCTION assignment in the output
  %     file. This does not affect the derivatives.
  
  % Check for yRow varying and yCol > 1
  if yOver.func.size(2) > 1 && sum(ySizes(1,:)~=ySizes(1,1))
    % the indices in ASGNINDS are no longer valid, need to convert
    % them to match the overmapped Y.
    NewInds = reshape(ASGNINDS,NUMinds,prod(FORLENGTHS));
    SScount = 0; MaxRef = zeros(yOver.func.size(1),yOver.func.size(2));
    MaxRef(:) = 1:numel(MaxRef);
    for Icount = 1:size(NewInds,2)
      if sum(NewInds(:,Icount))
        SScount = SScount+1;
        yRef = MaxRef(1:ySizes(1,SScount),1:ySizes(2,RScount));
        OldInd = NewInds(:,Icount);
        NewInds(logical(OldInd),Icount) = yRef(nonzeros(OldInd));
      end
    end
    ASGNINDS = reshape(NewInds,size(ASGNINDS,1),size(ASGNINDS,2));
    % ASGNINDS are now changed to match OverMap.
  end
  
  XSIZES = cell(2,1);
  XSIZEDEP = cell(2,1);
  % Check for changing x-size
  for SScount = 1:2
    xSize = xSizes(SScount,:);
    if sum(xSize ~= xSize(1))
      % This dimension changes - put it into same form as indice
      % array.
      xSizeTemp = zeros(1,prod(FORLENGTHS));
      xSizeTemp(logical(sum(reshape(ASGNINDS,NUMinds,...
        prod(FORLENGTHS)),1))) = xSize;
      xSize = reshape(xSizeTemp,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % Get SIZE dependencies
      XSIZEDEP{SScount} = GetDataDependencies(xSize,FORLENGTHS);
      % Remove any unneeded sizes
      XSIZES{SScount} = RemoveUnneededData(xSize,FORLENGTHS,...
        XSIZEDEP{SScount});
    end
  end
  
  %                    Get the Derivative Indices                         %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
  % This gets all of the derivative indices which need to get
  % printed. - GetSubsasgnInds does all the hard work.
  % Get SPARSE indice Dependencies
  [ASGNDEP] = GetDataDependencies(ASGNINDS,FORLENGTHS);
  % Remove any unneeded indices
  DEPINDS = RemoveUnneededData(ASGNINDS,FORLENGTHS,ASGNDEP);
  % Check for b changing size and bCols > 1
  if xOver.func.size(2) > 1 && sum(diff(xSizes(1,:)))
    % x changes row size and is a matrix, need to put xSizes into
    % same form as DEPINDS so that we can properly get the
    % derivative indices.
    xSizeTemp = zeros(2,prod(FORLENGTHS));
    xSizeTemp(:,logical(sum(reshape(ASGNINDS,NUMinds,prod(FORLENGTHS)),1))) = xSizes;
    xSizeTemp = reshape(xSizeTemp,2*prod(FORLENGTHS(2:end)),FORLENGTHS(1));
    xSizes = RemoveUnneededData(xSizeTemp,FORLENGTHS,ASGNDEP);
    if size(xSizes,2) > 1
      xSizes = reshape(xSizes,2,prod(FORLENGTHS(ASGNDEP)));
    end
  else
    xSizes = [];
  end
  % Get the Derivative indices
  if size(DEPINDS,2) > 1
    DEPINDS = reshape(DEPINDS,NUMinds,prod(FORLENGTHS(ASGNDEP)));
  end
  DerInds = GetSubsasgnInds(yOver,xOver,DEPINDS,xSizes);
  
  
  
  %                   Assign Names for Printing Run                       %
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

  % Pre-Allocate SPARSE INDICES cells - We store derivative index mapping
  % in .INDICES and input/outputs changing size in .SIZES
  for Fcount = 1:length(FORLENGTHS)
    ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SPARSE(ASGNLOCS(2,Fcount))...
      .INDICES = cell(NUMvod,3);
    ADIGATORFORDATA(ASGNLOCS(1,Fcount)).SPARSE(ASGNLOCS(2,Fcount))...
      .SIZES   = cell(2,3);
  end
  
  % Derivatives.
  for Vcount = 1:NUMvod
    if ~isempty(DerInds{Vcount,1})
      % These get stored in .INDICES{Vcount,1:3}
      DataNameCount = AssignReferenceNames(DerInds{Vcount,1},ASGNDEP,...
        ASGNLOCS,FORLENGTHS,'SPARSE','INDICES',Vcount,0,DataNameCount);
    end
  end
  
  % Function Sizes
  for SScount = 1:2
    if ~isempty(XSIZES{SScount})
      % These get stored in .SIZES{SScount,1:3}, where SScount is dimension
      % of input
      DataNameCount = AssignReferenceNames(XSIZES{SScount},XSIZEDEP{SScount},...
        ASGNLOCS,FORLENGTHS,'SPARSE','SIZES',SScount,0,DataNameCount);
    end
  end
end
% ----------------------------------------------------------------------- %
%                             NONZEROS                                    %
% ----------------------------------------------------------------------- %
for Ncount   = 1:length(ADIGATORFORDATA(FORCOUNT).NONZEROS)
  NZINDS     = ADIGATORFORDATA(FORCOUNT).NONZEROS(Ncount).INDICES;
  NZLOCS     = ADIGATORFORDATA(FORCOUNT).NONZEROS(Ncount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,NZLOCS);
  if isempty(NZINDS)
    % This call to NONZEROS never fires
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(NZLOCS(1,Fcount)).NONZEROS(NZLOCS(2,Fcount))...
        .INDICES = cell(1,3);
      ADIGATORFORDATA(NZLOCS(1,Fcount)).NONZEROS(NZLOCS(2,Fcount))...
        .SIZES   = [];
    end
    continue
  elseif size(NZINDS,2) < OuterLoopMaxLength
    NZINDS(end,OuterLoopMaxLength) = 0;
  end
  NZFLAG = ADIGATORFORDATA(NZLOCS(1,end)).NONZEROS(NZLOCS(2,end)).FLAGS(1);

  % Maximum number of references per NONZEROS call
  NUMinds = numel(NZINDS)/prod(FORLENGTHS);
  ADIGATORFORDATA(FORCOUNT).NONZEROS(Ncount).FLAGS(2) = NUMinds;
  
  nzIndPrint = 0;
  % Two Cases where we need to print out Function indices:
  %     1. When x is a matrix and its row size changes
  %     2. When x has non-zero entries.
  
  % Check to see if x is a matrix and row size changes
  xSizes = ADIGATORFORDATA(NZLOCS(1,end)).NONZEROS(NZLOCS(2,end)).SIZES;
  if max(xSizes(2,:)) > 1 && sum(xSizes(1,:)~=xSizes(1,1))
    % X is matrix, row size changes, need to print out indices for
    % every case.
    % Re-Write indices so they are proper to overmapped x
    xRowMax = max(xSizes(1,:)); xColMax = max(xSizes(2,:));
    xRefMax = zeros(xRowMax,xColMax); xRefMax(:) = 1:xRowMax*xColMax;
    NewInds = reshape(NZINDS,NUMinds,prod(FORLENGTHS));
    NScount = 0;
    for NIcount = 1:size(NewInds,2)
      if sum(NewInds(:,NIcount))
        NScount = NScount+1;
        xRef = xRefMax(1:xSizes(1,NScount),1:xSizes(2,NScount));
        OldInds = NewInds(:,NIcount);
        NewInds(logical(OldInds),NIcount) = xRef(nonzeros(OldInds));
      end
    end
    NZINDS = reshape(NewInds,size(NZINDS,1),size(NZINDS,2));
    nzIndPrint = 1;
  elseif NZFLAG == 3
    % x is symbolic and has non-zero entries, need to do
    % referencing
    nzIndPrint = 1;
  end
  
  % Pre-Allocate NONZERO INDICES cells - here we use .INDICES field to
  % store function reference indices - clear .SIZES field
  for Fcount = 1:length(FORLENGTHS)
    ADIGATORFORDATA(NZLOCS(1,Fcount)).NONZEROS(NZLOCS(2,Fcount))...
      .INDICES = cell(1,3);
    ADIGATORFORDATA(NZLOCS(1,Fcount)).NONZEROS(NZLOCS(2,Fcount))...
      .SIZES   = [];
  end
  
  if nzIndPrint
    % Get NONZEROS indice Dependencies
    [NZDEP] = GetDataDependencies(NZINDS,FORLENGTHS);
    % Remove any unneeded indices
    DEPINDS = RemoveUnneededData(NZINDS,FORLENGTHS,NZDEP);
    
    % Assign Indice names that will get printed
    DataNameCount = AssignReferenceNames(DEPINDS,NZDEP,NZLOCS,FORLENGTHS,...
      'NONZEROS','INDICES',1,0,DataNameCount);
  end
end
% ----------------------------------------------------------------------- %
%                              HORZCAT                                    %
% ----------------------------------------------------------------------- %
for Hcount   = 1:length(ADIGATORFORDATA(FORCOUNT).HORZCAT)
  % y = horzcat(Inputs);
  HSIZES     = ADIGATORFORDATA(FORCOUNT).HORZCAT(Hcount).SIZES;
  HLOCS      = ADIGATORFORDATA(FORCOUNT).HORZCAT(Hcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,HLOCS);
  if isempty(HSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .INDICES = [];
    end
    continue
  elseif size(HSIZES,2) < OuterLoopMaxLength
    HSIZES(end,OuterLoopMaxLength) = 0;
  end
  
  % Grab the OverMapped variables
  HVARS      = ADIGATORFORDATA(HLOCS(1,end)).HORZCAT(HLOCS(2,end)).VARS;
  NUMHinput  = size(HVARS,1)-1;
  for Icount = 2:NUMHinput+1
    if ~isa(HVARS{Icount},'cada')
      HVARS{Icount} = ADIGATORVARIABLESTORAGE.OVERMAP{HVARS{Icount}};
    end
  end
  ADIGATORFORDATA(HLOCS(1,end)).HORZCAT(HLOCS(2,end)).VARS = HVARS;
  
  % ---See if Sizes of Y or Inputs Change---%
  [HORZDEP] = GetDataDependencies(HSIZES,FORLENGTHS);

  if sum(HORZDEP) % The sizes do change somewhere
    % ---Pre-Allocate HORZCAT Size and Indices Cells - we use .INDICES to 
    % store derivative indices and .SIZES to store input/output sizes --- %
    for Fcount = 1:length(HORZDEP)
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .SIZES   = cell(NUMHinput,3);
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,NUMHinput,3);
    end
    
    % Size Changes Somewhere - Remove the Independent Sizes
    [HDEPSIZES] = RemoveUnneededData(HSIZES,FORLENGTHS,HORZDEP);
    
    if size(HDEPSIZES,2) > 1
      HDEPSIZES = reshape(HDEPSIZES,1+NUMHinput,prod(FORLENGTHS(HORZDEP)));
    end
    
    % ------------Get the Derivative Indices--------------------- %
    HDerInds = GetHorzcatInds(HDEPSIZES,HVARS);

    % Assign Derivative Indice names that will get printed
    for Vcount = 1:NUMvod
      for Icount = 1:NUMHinput
        DerInds = HDerInds{Vcount,Icount};
        % These Derivative Indices get assigned to
        % .INDICES{Vcount,Icount,1:3}
        if ~isempty(DerInds)
          DataNameCount = AssignReferenceNamesCat(DerInds,HORZDEP,HLOCS,...
            FORLENGTHS,'HORZCAT','INDICES',[Vcount Icount],DataNameCount);
        end
      end
    end
    
    % -------------Function Sizes Changing------------------------- %
    HSIZES = reshape(HSIZES,NUMHinput+1,prod(FORLENGTHS));
    for Icount = 1:NUMHinput
      ISIZES = HSIZES(Icount+1,:);
      ISIZES = reshape(ISIZES,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % These dont all necessarily have the same dependence as the
      % dependence we use for the derivatives, so need to do each one
      % seperately.
      HORZDEP = GetDataDependencies(ISIZES,FORLENGTHS);
      
      if sum(HORZDEP)
        IndLengths = RemoveUnneededData(ISIZES,FORLENGTHS,HORZDEP);
        DataNameCount = AssignReferenceNames(IndLengths,HORZDEP,HLOCS,...
          FORLENGTHS,'HORZCAT','SIZES',Icount,0,DataNameCount);
      end
    end
  else
    % The sizes dont change - Clear .SIZES and .INDICES
    for Fcount = 1:length(HORZDEP)
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(HLOCS(1,Fcount)).HORZCAT(HLOCS(2,Fcount))...
        .INDICES = [];
    end
  end
end
% ----------------------------------------------------------------------- %
%                              VERTCAT                                    %
% ----------------------------------------------------------------------- %
for VEcount  = 1:length(ADIGATORFORDATA(FORCOUNT).VERTCAT)
  % y = horzcat(Inputs);
  VSIZES     = ADIGATORFORDATA(FORCOUNT).VERTCAT(VEcount).SIZES;
  VLOCS      = ADIGATORFORDATA(FORCOUNT).VERTCAT(VEcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,VLOCS);
  if isempty(VSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .INDICES = [];
    end
    continue
  elseif size(VSIZES,2) < OuterLoopMaxLength
    VSIZES(end,OuterLoopMaxLength) = 0;
  end
  % Grab the OverMapped variables
  VVARS     = ADIGATORFORDATA(VLOCS(1,end)).VERTCAT(VLOCS(2,end)).VARS;
  NUMVinput = size(VVARS,1)-1;
  for Icount = 2:NUMVinput+1
    if ~isa(VVARS{Icount},'cada')
      VVARS{Icount} = ADIGATORVARIABLESTORAGE.OVERMAP{VVARS{Icount}};
    end
  end
  ADIGATORFORDATA(VLOCS(1,end)).VERTCAT(VLOCS(2,end)).VARS = VVARS;
  
  % ---See if Sizes of Y or Inputs Change---%
  [VERTDEP] = GetDataDependencies(VSIZES,FORLENGTHS);

  if sum(VERTDEP) % Sizes change somewhere
    % ---Pre-Allocate HORZCAT Size and Indices Cells - we use .INDICES to 
    % store derivative indices and .SIZES to store input/output sizes --- %
    for Fcount = 1:length(VERTDEP)
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .SIZES   = cell(NUMVinput,3);
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,NUMVinput,3);
    end
    % Remove the Independent Sizes
    [VDEPSIZES] = RemoveUnneededData(VSIZES,FORLENGTHS,VERTDEP);
    if size(VDEPSIZES,2) > 1
      VDEPSIZES = reshape(VDEPSIZES,1+NUMVinput,prod(FORLENGTHS(VERTDEP)));
    end
    
    % ------------Get the Derivative Indices--------------------- %
    VDerInds = GetVertcatInds(VDEPSIZES,VVARS);
    % Assign Derivative Indice names that will get printed
    for Vcount = 1:NUMvod
      for Icount = 1:NUMVinput
        DerInds = VDerInds{Vcount,Icount};
        if ~isempty(DerInds)
          % These Derivative Indices get assigned to
          % .INDICES{Vcount,Icount,1:3}
          DataNameCount = AssignReferenceNamesCat(DerInds,VERTDEP,VLOCS,...
            FORLENGTHS,'VERTCAT','INDICES',[Vcount Icount],DataNameCount);
        end
      end
    end
    
    % -------------Function Sizes Changing------------------------- %
    VSIZES = reshape(VSIZES,NUMVinput+1,prod(FORLENGTHS));
    for Icount = 1:NUMVinput
      ISIZES = VSIZES(Icount+1,:);
      ISIZES = reshape(ISIZES,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % These dont all necessarily have the same dependence as the
      % dependence we use for the derivatives, so need to do each one
      % seperately.
      VERTDEP = GetDataDependencies(ISIZES,FORLENGTHS);
      if sum(VERTDEP)
        IndLengths = RemoveUnneededData(ISIZES,FORLENGTHS,VERTDEP);
        DataNameCount = AssignReferenceNames(IndLengths,VERTDEP,VLOCS,...
          FORLENGTHS,'VERTCAT','SIZES',Icount,0,DataNameCount);
      end
    end
  else
    % The sizes dont change - Clear .SIZES and .INDICES
    for Fcount = 1:length(VERTDEP)
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(VLOCS(1,Fcount)).VERTCAT(VLOCS(2,Fcount))...
        .INDICES = [];
    end
  end
end
% ----------------------------------------------------------------------- %
%                             TRANSPOSE                                   %
% ----------------------------------------------------------------------- %
for Tcount   = 1:length(ADIGATORFORDATA(FORCOUNT).TRANSPOSE)
  % y = transpose(x)
  TSIZES     = ADIGATORFORDATA(FORCOUNT).TRANSPOSE(Tcount).SIZES;
  TRANLOCS   = ADIGATORFORDATA(FORCOUNT).TRANSPOSE(Tcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,TRANLOCS);
  if isempty(TSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .INDICES = [];
    end
    continue
  elseif size(TSIZES,2) < OuterLoopMaxLength
    TSIZES(end,OuterLoopMaxLength) = 0;
  end
  % Grab the OverMapped variables
  yOver = ADIGATORFORDATA(TRANLOCS(1,end)).TRANSPOSE(TRANLOCS(2,end)).VARS{1};
  xOver = ADIGATORFORDATA(TRANLOCS(1,end)).TRANSPOSE(TRANLOCS(2,end)).VARS{2};
  if ~isa(xOver,'cada')
    xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
    ADIGATORFORDATA(TRANLOCS(1,end)).TRANSPOSE(TRANLOCS(2,end)).VARS{2} = xOver;
  end
  
  % ---See if Sizes of Y or X Change---%
  [TRANDEP] = GetDataDependencies(TSIZES,FORLENGTHS);
  
  if sum(TRANDEP) && xOver.func.size(1) > 1 && xOver.func.size(2) > 1 && ...
      ~any(isinf(xOver.func.size))
    % ---Pre-Allocate TRANSPOSE Size and Indice Cells--- %
    for Fcount = 1:length(TRANDEP)
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .SIZES   = cell(2,3);
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,3);
    end
    
    % Size Changes Somewhere - Remove the Independent Sizes
    [TDEPSIZES] = RemoveUnneededData(TSIZES,FORLENGTHS,TRANDEP);
    if size(TDEPSIZES,2) > 1
      TDEPSIZES = reshape(TDEPSIZES,2,prod(FORLENGTHS(TRANDEP)));
    end
    
    % ------------Get the Derivative Indices--------------------- %
    TDerInds = GetTransposeInds(TDEPSIZES,yOver,xOver);
    
    % Assign Derivative Indice names that will get printed
    for Vcount = 1:NUMvod
      DerInds = TDerInds{Vcount};
      if ~isempty(DerInds)
        DataNameCount = AssignReferenceNames(DerInds,TRANDEP,TRANLOCS,...
          FORLENGTHS,'TRANSPOSE','INDICES',Vcount,0,DataNameCount);
      end
    end
    
    % -------------Function Sizes Changing----------------------- %
    TSIZES = reshape(TSIZES,2,prod(FORLENGTHS));
    for Icount = 1:2
      ISIZES = TSIZES(Icount,:);
      ISIZES = reshape(ISIZES,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % These dont all necessarily have the same dependence as the
      % dependence we use for the derivatives, so need to do each one
      % seperately.
      TRANDEP = GetDataDependencies(ISIZES,FORLENGTHS);
      
      if sum(TRANDEP)
        IndLengths    = RemoveUnneededData(ISIZES,FORLENGTHS,TRANDEP);
        DataNameCount = AssignReferenceNames(IndLengths,TRANDEP,TRANLOCS,...
          FORLENGTHS,'TRANSPOSE','SIZES',Icount,0,DataNameCount);
      end
    end
  else
    % No size change - Clear .SIZES and .INDICES
    for Fcount = 1:length(TRANDEP)
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(TRANLOCS(1,Fcount)).TRANSPOSE(TRANLOCS(2,Fcount))...
        .INDICES = [];
    end
  end
end
% ----------------------------------------------------------------------- %
%                               REPMAT                                    %
% ----------------------------------------------------------------------- %
for REPcount = 1:length(ADIGATORFORDATA(FORCOUNT).REPMAT)
  % y = repmat(x,M,N);
  RSIZES     = ADIGATORFORDATA(FORCOUNT).REPMAT(REPcount).SIZES;
  REPLOCS    = ADIGATORFORDATA(FORCOUNT).REPMAT(REPcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,REPLOCS);
  if isempty(RSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .INDICES = [];
    end
    continue
  elseif size(RSIZES,2) < OuterLoopMaxLength
    RSIZES(end,OuterLoopMaxLength) = 0;
  end
  % Grab the OverMapped variables
  yOver = ADIGATORFORDATA(REPLOCS(1,end)).REPMAT(REPLOCS(2,end)).VARS{1};
  xOver = ADIGATORFORDATA(REPLOCS(1,end)).REPMAT(REPLOCS(2,end)).VARS{2};
  if ~isa(xOver,'cada')
    xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
    ADIGATORFORDATA(REPLOCS(1,end)).REPMAT(REPLOCS(2,end)).VARS{2} = xOver;
  end
  
  % ---See if Sizes of Y or X Change---%
  [REPDEP] = GetDataDependencies(RSIZES,FORLENGTHS);
  
  if sum(REPDEP)
    % ---Pre-Allocate REPMAT Size and Indice Cells--- %
    for Fcount = 1:length(REPDEP)
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .SIZES   = cell(2,3);
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,3);
    end
    
    % Size Changes Somewhere - Remove the Independent Sizes
    [RDEPSIZES] = RemoveUnneededData(RSIZES,FORLENGTHS,REPDEP);
    if size(RDEPSIZES,2) > 1
      RDEPSIZES = reshape(RDEPSIZES,4,prod(FORLENGTHS(REPDEP)));
    end
    
    % ------------Get the Derivative Indices--------------------- %
    RDerInds = GetRepmatInds(RDEPSIZES,yOver,xOver);
    
    % Assign Derivative Indice names that will get printed
    for Vcount = 1:NUMvod
      DerInds = RDerInds{Vcount};
      if ~isempty(DerInds)
        DataNameCount = AssignReferenceNames(DerInds,REPDEP,REPLOCS,...
          FORLENGTHS,'REPMAT','INDICES',Vcount,0,DataNameCount);
      end
    end
    
    % -------------Function Sizes Changing----------------------- %
    RSIZES = reshape(RSIZES,4,prod(FORLENGTHS));
    RSIZES = RSIZES(3:4,:);
    for Icount = 1:2
      ISIZES = RSIZES(Icount,:);
      ISIZES = reshape(ISIZES,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % These dont all necessarily have the same dependence as the
      % dependence we use for the derivatives, so need to do each one
      % seperately.
      REPDEP = GetDataDependencies(ISIZES,FORLENGTHS);
      
      if sum(REPDEP)
        IndLengths = RemoveUnneededData(ISIZES,FORLENGTHS,REPDEP);
        DataNameCount = AssignReferenceNames(IndLengths,REPDEP,REPLOCS,...
          FORLENGTHS,'REPMAT','SIZES',Icount,0,DataNameCount);
      end
    end
  else
    % Clear SIZES and INDICES fields
    for Fcount = 1:length(REPDEP)
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(REPLOCS(1,Fcount)).REPMAT(REPLOCS(2,Fcount))...
        .INDICES = [];
    end
  end
end
% ----------------------------------------------------------------------- %
%                              RESHAPE                                    %
% ----------------------------------------------------------------------- %
for REScount = 1:length(ADIGATORFORDATA(FORCOUNT).RESHAPE)
  % y = reshape(x,M,N);
  RSIZES     = ADIGATORFORDATA(FORCOUNT).RESHAPE(REScount).SIZES;
  RESLOCS    = ADIGATORFORDATA(FORCOUNT).RESHAPE(REScount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,RESLOCS);
  if isempty(RSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .INDICES = [];
    end
    continue
  elseif size(RSIZES,2) < OuterLoopMaxLength
    RSIZES(end,OuterLoopMaxLength) = 0;
  end
  % Grab the OverMapped variables
  yOver = ADIGATORFORDATA(RESLOCS(1,end)).RESHAPE(RESLOCS(2,end)).VARS{1};
  xOver = ADIGATORFORDATA(RESLOCS(1,end)).RESHAPE(RESLOCS(2,end)).VARS{2};
  if ~isa(xOver,'cada')
    xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
    ADIGATORFORDATA(RESLOCS(1,end)).RESHAPE(RESLOCS(2,end)).VARS{2} = xOver;
  end
  
  % ---See if Sizes of Y or X Change---%
  [RESDEP] = GetDataDependencies(RSIZES,FORLENGTHS);
  
  if sum(RESDEP)
    % ---Pre-Allocate RESHAPE Size and Indice Cells--- %
    for Fcount = 1:length(RESDEP)
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .SIZES   = cell(2,3);
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .INDICES = cell(NUMvod,3);
    end
    
    % Size Changes Somewhere - Remove the Independent Sizes
    [RDEPSIZES] = RemoveUnneededData(RSIZES,FORLENGTHS,RESDEP);
    if size(RDEPSIZES,2) > 1
      RDEPSIZES = reshape(RDEPSIZES,4,prod(FORLENGTHS(RESDEP)));
    end
    
    % ------------Get the Derivative Indices--------------------- %
    RDerInds = GetReshapeInds(RDEPSIZES,yOver,xOver);
    
    % Assign Derivative Indice names that will get printed
    for Vcount = 1:NUMvod
      DerInds = RDerInds{Vcount};
      if ~isempty(DerInds)
        DataNameCount = AssignReferenceNames(DerInds,RESDEP,RESLOCS,...
          FORLENGTHS,'RESHAPE','INDICES',Vcount,0,DataNameCount);
      end
    end
    
    % -------------Function Sizes Changing----------------------- %
    RSIZES = reshape(RSIZES,4,prod(FORLENGTHS));
    RSIZES = RSIZES(3:4,:);
    for Icount = 1:2
      ISIZES = RSIZES(Icount,:);
      ISIZES = reshape(ISIZES,prod(FORLENGTHS(2:end)),FORLENGTHS(1));
      % These dont all necessarily have the same dependence as the
      % dependence we use for the derivatives, so need to do each one
      % seperately.
      RESDEP = GetDataDependencies(ISIZES,FORLENGTHS);
      
      if sum(RESDEP)
        IndLengths = RemoveUnneededData(ISIZES,FORLENGTHS,RESDEP);
        DataNameCount = AssignReferenceNames(IndLengths,RESDEP,RESLOCS,...
          FORLENGTHS,'RESHAPE','SIZES',Icount,0,DataNameCount);
      end
    end
  else
    % Clear SIZES and INDICES fields
    for Fcount = 1:length(RESDEP)
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .SIZES   = [];
      ADIGATORFORDATA(RESLOCS(1,Fcount)).RESHAPE(RESLOCS(2,Fcount))...
        .INDICES = [];
    end
  end
end
% ----------------------------------------------------------------------- %
%                               SIZE                                      %
% ----------------------------------------------------------------------- %
for SIZcount = 1:length(ADIGATORFORDATA(FORCOUNT).SIZE)
  % This has data from SIZE,LENGTH,NUMEL - just need to print out
  % their outputs if they are changing.
  OUTSIZES   = ADIGATORFORDATA(FORCOUNT).SIZE(SIZcount).SIZES;
  SIZLOCS    = ADIGATORFORDATA(FORCOUNT).SIZE(SIZcount).LOCS;
  FORLENGTHS = GetForLengths(ADIGATORFORDATA,SIZLOCS);
  if isempty(OUTSIZES)
    for Fcount = 1:length(FORLENGTHS)
      ADIGATORFORDATA(SIZLOCS(1,Fcount)).SIZE(SIZLOCS(2,Fcount))...
        .SIZES = [];
    end
    continue
  elseif size(OUTSIZES,2) < OuterLoopMaxLength
    OUTSIZES(end,OuterLoopMaxLength) = 0;
  end
  OUTSIZES(isinf(OUTSIZES)) = -1;
  % ---See if Output Sizes Change---%
  [SIZDEP]   = GetDataDependencies(OUTSIZES,FORLENGTHS);
  OUTSIZES(OUTSIZES == -1) = Inf;
  OUTSIZES(isnan(OUTSIZES)) = 0;
  
  [DEPSIZES] = RemoveUnneededData(OUTSIZES,FORLENGTHS,SIZDEP);
  if sum(SIZDEP)
    % ---Pre-Allocate SIZE Size and Indice Cells--- %
    for Fcount = 1:length(SIZDEP)
      ADIGATORFORDATA(SIZLOCS(1,Fcount)).SIZE(SIZLOCS(2,Fcount))...
        .SIZES = cell(1,3);
    end
    
    % Size Changes Somewhere - Remove the Independent Sizes
    [DEPSIZES] = RemoveUnneededData(OUTSIZES,FORLENGTHS,SIZDEP);
    NUMoutput = numel(DEPSIZES)/prod(FORLENGTHS(SIZDEP));
    if size(DEPSIZES,2) > 1
      DEPSIZES = reshape(DEPSIZES,NUMoutput,prod(FORLENGTHS(SIZDEP)));
    end
    
    % -------------Function Sizes Changing----------------------- %
    DataNameCount = AssignReferenceNames(DEPSIZES,SIZDEP,SIZLOCS,...
      FORLENGTHS,'SIZE','SIZES',1,0,DataNameCount);
  else
    % clear SIZES field - store the DEPSIZES
    for Fcount = 1:length(SIZDEP)-1
      ADIGATORFORDATA(SIZLOCS(1,Fcount)).SIZE(SIZLOCS(2,Fcount))...
        .SIZES = [];
    end
    ADIGATORFORDATA(SIZLOCS(1,end)).SIZE(SIZLOCS(2,end)).SIZES = DEPSIZES;
  end
end

return
end

function FORLENGTHS = GetForLengths(ADIGATORFORDATA,FunLocs)
NUMloops   = size(FunLocs,2);
FORLENGTHS = zeros(1,NUMloops);
for Fcount = 1:NUMloops
  FORLENGTHS(Fcount) = ADIGATORFORDATA(FunLocs(1,Fcount)).MAXLENGTH;
end

end

function Dependency = GetDataDependencies(Inds,Lengths)
% Function recursively calls itself to see what indices are dependent upon
% what FOR loops
NUMloops = length(Lengths);
Dependency = zeros(1,NUMloops);
% OuterLoop Dependency

if diff(Inds,1,2) == 0
  Dependency(1) = 0;
else
  NUMinds = size(Inds,2);
  ZeroMap = ones(NUMinds,1);
  for Icount = 1:NUMinds
    if Inds(:,Icount) == 0
      ZeroMap(Icount) = 0;
    end
  end
  NZinds = Inds(:,logical(ZeroMap));
  if size(NZinds,2) == 1
    Dependency(1) = 0;
  elseif diff(NZinds,1,2)==0
    Dependency(1) = 0;
  else
    Dependency(1) = 1;
  end
end

% Children loops
if NUMloops > 1
  Clength = size(Inds,1)/Lengths(2);
  CLengths = Lengths(2:end);
  for Icount = 1:Lengths(1)
    Cinds = Inds(:,Icount);
    Cinds = reshape(Cinds,Clength,Lengths(2));
    CDep = GetDataDependencies(Cinds,CLengths);
    Dependency(2:end) = Dependency(2:end)+CDep;
  end
end
Dependency = logical(Dependency);
return
end

function NewInds = RemoveUnneededData(Inds,Lengths,Dependency)
NUMloops = length(Lengths);
% Indices need removed
% Get to the first dependent loop
if NUMloops == 1
  if ~Dependency(1)
    NewInds = [];
    for Icount = 1:size(Inds,2)
      if any(Inds(:,Icount))
        NewInds = Inds(:,Icount);
        break
      end
    end
    if isempty(NewInds)
      NewInds = Inds(:,1);
    end
  else
    NewInds = Inds;
  end
  return
end

Dcount = 1;
while ~Dependency(Dcount)
  Dcount = Dcount+1;
  for Icount = 1:size(Inds,2)
    if any(Inds(:,Icount))
      Inds = Inds(:,Icount);
      break
    end
  end
  Inds = reshape(Inds,length(Inds)/Lengths(Dcount),Lengths(Dcount));
  if Dcount == NUMloops
    if ~Dependency(Dcount)
      for Icount = 1:size(Inds,2)
        if any(Inds(:,Icount))
          Inds = Inds(:,Icount);
          break
        end
      end
    end
    break
  end
end

if Dcount < NUMloops && sum(Dependency(Dcount:end)) < length(Dependency(Dcount:end))
  % is not the Innermost loop, and the inner loops have an independency
  % somewhere
  NewInds = zeros(1,Lengths(Dcount));
  IndLength = size(Inds,1);
  IndLength = IndLength/Lengths(Dcount+1);
  Inds1 = reshape(Inds(:,1),IndLength,Lengths(Dcount+1));
  Inds1 = RemoveUnneededData(Inds1,Lengths(Dcount+1:end),Dependency(Dcount+1:end));
  NewInds(1:numel(Inds1),1) = Inds1(:);
  for Icount = 2:Lengths(Dcount)
    Inds1 = reshape(Inds(:,Icount),IndLength,Lengths(Dcount+1));
    Inds1 = RemoveUnneededData(Inds1,Lengths(Dcount+1:end),Dependency(Dcount+1:end));
    NewInds(:,Icount) = Inds1(:);
  end
else
  NewInds = Inds;
end
return
end

function DerInds = GetSubsrefInds(x,y,inds,ySizes)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;

xMrow = x.func.size(1); xNcol = x.func.size(2);
yMrow = y.func.size(1); yNcol = y.func.size(2);
if isinf(xMrow); xMrow = 1; elseif isinf(xNcol); xNcol = 1; end
if isinf(yMrow); yMrow = 1; elseif isinf(yNcol); yNcol = 1; end
yRef = zeros(yMrow,yNcol); yRef(:) = 1:yMrow*yNcol;

[~,DPFLAG] = cadafuncname(y.id);
ForLength = size(inds,2);
DerInds = cell(NUMvod,1);
if DPFLAG
  % Get the Derivative Indices - is y = x(inds), x and y OverMapped
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs) && ~isempty(x.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      % get dx deriv pattern
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2);
      nzx   = length(xrows);
      % dx map with ones
      dx = sparse(xrows,xcols,2*ones(nzx,1),xMrow*xNcol,nv);
      % dx map with index
      dxi = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
      % get dy deriv pattern
      yrows = y.deriv(Vcount).nzlocs(:,1);
      ycols = y.deriv(Vcount).nzlocs(:,2);
      nzy = length(yrows);
      % dy map with ones
      dyover = sparse(yrows,ycols,ones(nzy,1),yMrow*yNcol,nv);
      % dy map with index
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % Pre-Allocate Derivative Indices
      DerInds{Vcount} = zeros(nzy,ForLength);
      for Icount = 1:ForLength
        isubs = nonzeros(inds(:,Icount));
        % Number of indices may change - get the dy and dyi for this
        % particular length of indices
        NUMsubs = length(isubs);
        if ~isempty(ySizes)
          yiref = yRef(1:ySizes(1,Icount),1:ySizes(2,Icount));
          ysub = yiref(1:NUMsubs);
        else
          ysub = 1:NUMsubs;
        end
        dy = dyover(ysub,:);
        dyi = dyoveri(ysub,:);
        yind = nonzeros(dyi);
        
        % We want to find the locations in dy and the corresponding index
        % locations from dx that overlap
        dxlocsi = nonzeros(dxi(isubs,:));
        dxlocs = nonzeros(dy+dx(isubs,:));
        yindlocs = dxlocs(dxlocs==3 | dxlocs==1);
        yindlocs = yind(yindlocs == 3);
        
        dxlocs = dxlocs(dxlocs>1);
        xindlocs = dxlocsi(dxlocs == 3);
        
        DerInds{Vcount}(yindlocs,Icount) = xindlocs;
      end
    end
  end
end

return
end

function DerInds = GetSubsasgnInds(x,b,inds,bSizes)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;

xMrow = x.func.size(1); xNcol = x.func.size(2);
bMrow = b.func.size(1); bNcol = b.func.size(2);
if isinf(xMrow); xMrow = 1; elseif isinf(xNcol); xNcol = 1; end
if isinf(bMrow); bMrow = 1; elseif isinf(bNcol); bNcol = 1; end
bRef = zeros(bMrow,bNcol); bRef(:) = 1:bMrow*bNcol;

[~,DPFLAG] = cadafuncname(x.id);

ForLength = size(inds,2);
DerInds = cell(NUMvod,2);
if DPFLAG
  for Vcount = 1:ADIGATOR.NVAROFDIFF
    % if b has derivatives, x has to (overmapped)
    if ~isempty(x.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2);
      nzx = length(xrows);
      % dx map with ones
      dx = sparse(xrows,xcols,2*ones(nzx,1),xMrow*xNcol,nv);
      % dx map with index
      dxi = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
      if ~isempty(b.deriv(Vcount).nzlocs)
        brows = b.deriv(Vcount).nzlocs(:,1);
        bcols = b.deriv(Vcount).nzlocs(:,2);
        % db overmap with index
        dboveri = sparse(brows,bcols,1:length(brows),bMrow*bNcol,nv);
        DerInds{Vcount,1} = zeros(nzx,ForLength);
        DerInds{Vcount,2} = zeros(nzx,ForLength);
        for Icount = 1:ForLength
          isubs = nonzeros(inds(:,Icount));
          NUMsubs = length(isubs);
          if ~isempty(bSizes)
            biref = bRef(1:bSizes(1,Icount),1:bSizes(2,Icount));
            bsubs = biref(1:NUMsubs);
          else
            bsubs = 1:NUMsubs;
          end
          % db for this particular iteration with index
          dbi = dboveri(bsubs,:);
          
          [brows,bcols,bind] = find(dbi);
          % db for this particular iteration with ones
          db = sparse(brows,bcols,ones(size(brows)),NUMsubs,nv);
          
          dxlocsi = nonzeros(dxi(isubs,:));
          dxlocs = nonzeros(db + dx(isubs,:));
          bindlocs = dxlocs(dxlocs==1 | dxlocs==3);
          bind = bind(bindlocs == 3);
          dxlocs = dxlocs(dxlocs>1);
          dblocs = dxlocsi(dxlocs==3);
          zlocs = dxlocsi(dxlocs==2);
          % dblocs are the locations of the derivatives which are
          % assigned on this iteration
          % bind are the locations of the derivatives of the
          % b derivatives which are assigned on this iteration
          % zlocs are the locations of the derivatives which are
          % removed on this iteration
          DerInds{Vcount,1}(dblocs,Icount) = bind;
          DerInds{Vcount,2}(zlocs,Icount) = 1;
        end
      else
        DerInds{Vcount,2} = zeros(nzx,ForLength);
        for Icount = 1:ForLength
          isubs = nonzeros(inds(:,Icount));
          dxlocsi = nonzeros(dxi(isubs,:));
          DerInds{Vcount,2}(dxlocsi,Icount) = 1;
        end
      end
    end
  end
end
return
end

function DerInds = GetHorzcatInds(Sizes,Vars)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
NUMinput = length(Vars)-1;
NUMloops = size(Sizes,2);
y = Vars{1};
yMrow = y.func.size(1); yNcol = y.func.size(2);
if isinf(yMrow); yMrow = 1; end
yOverRef = zeros(yMrow,yNcol); yOverRef(:) = 1:yMrow*yNcol;
Inputs = Vars(2:end);

DerInds = cell(NUMvod,NUMinput);
[~,DPflag] = cadafuncname(y.id);

if DPflag
  
  yRowSizes = Sizes(1,:);
  yColSizes = sum(Sizes(2:end,:),1);
  iColSizes = Sizes(2:end,:);
  iColSums = zeros(NUMinput,NUMloops);
  for Icount = 1:NUMinput
    iColSums(Icount,:) = sum(iColSizes(1:Icount,:),1);
  end
  
  % See if Y is changing Row Size
  yRowsTemp = nonzeros(yRowSizes);
  if sum(yRowsTemp ~= yRowsTemp(1)); yRowDep = 1; else yRowDep = 0; end
  % If y is changing row size then everything is changing size
  
  if ~yRowDep
    % if y rows are changing, then check to see what is changing.
    yColsTemp = nonzeros(yColSizes);
    if sum(yColsTemp ~= yColsTemp(1)); yColDep = 1; else yColDep = 0; end
    
    iColsDep = zeros(NUMinput,1);
    for Icount = 1:NUMinput
      iColsTemp = nonzeros(Sizes(Icount+1,:));
      if sum(iColsTemp ~= iColsTemp(1)); iColsDep(Icount) = 1; end
    end
  end
  
  
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      yrows = y.deriv(Vcount).nzlocs(:,1); ycols = y.deriv(Vcount).nzlocs(:,2);
      nzy = length(yrows);
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % Overmapped DY with linear index of DY as the values.
      for Icount = 1:NUMinput
        x = Inputs{Icount};
        % X = Input(I);
        if ~isempty(x.deriv(Vcount).nzlocs);
          DerInds{Vcount,Icount} = zeros(nzy,NUMloops);
          xMrow = x.func.size(1); xNcol = x.func.size(2);
          if isinf(xMrow); xMrow = 1; end
          xrows = x.deriv(Vcount).nzlocs(:,1); xcols = x.deriv(Vcount).nzlocs(:,2);
          nzx = length(xrows);
          dxoveri = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
          % Overmapped DX with linear index of DX as the values.
          xOverRef = zeros(xMrow,xNcol); xOverRef(:) = 1:xMrow*xNcol;
          for Fcount = 1:NUMloops
            if yRowSizes(Fcount)
              % This loop has a horzcat in it.
              if yRowDep || (yColDep && iColsDep(Icount))
                % Both x and y are changing size
                yref = yOverRef(1:yRowSizes(Fcount),1:yColSizes(Fcount));
                % YREF - Function indices corresponding to Y on this loop.
                try
                xref = xOverRef(1:yRowSizes(Fcount),1:iColSizes(Icount,Fcount));
                catch
                  keyboard
                end
                % XREF - Function indices corresponding to X on this loop.
              elseif yColDep
                % Y is changing size, this input is not.
                xref = xOverRef;
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef(1:yRowSizes(Fcount),1:yColSizes(Fcount));
                % YREF - Function indices corresponding to Y on this loop.
              elseif iColsDep(Icount)
                % Y is not changing size, but this input is.
                xref = xOverRef(1:yRowSizes(Fcount),1:iColSizes(Icount,Fcount));
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef;
                % YREF - Function indices corresponding to Y on this loop.
              else
                % Neither Y nor X is changing size
                xref = xOverRef;
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef;
                % YREF - Function indices corresponding to Y on this loop.
              end
              dxi = dxoveri(xref(:),:);
              % DXI - Section of Overmapped DX corresponding to this loop -
              % with the linear index of DX as the values
              
              % Need to look at the part of DYI that DX is getting mapped
              % into
              if Icount == 1
                XintoYref = yref(:,1:iColSizes(1,Fcount));
              else
                XintoYref = yref(:,iColSums(Icount-1,Fcount)+1:iColSums(Icount,Fcount));
              end
              % XintoYref - Function indices of Overmapped Y which the X in
              % this loop is getting mapped into.
              dyi = dyoveri(XintoYref(:),:);
              % DYI - Section of Overmapped DY which corresponds to where
              % DXI is getting mapped into in this loop.
              
              % Now have DXI and DYI, need to see what derivatives are
              % common to both of them.
              [yrows,ycols] = find(dyi); [xrows,xcols] = find(dxi);
              dy1 = sparse(yrows,ycols,ones(size(yrows)),size(dyi,1),size(dyi,2));
              dx1 = sparse(xrows,xcols,ones(size(xrows)),size(dxi,1),size(dxi,2));
              dxy2 = dy1+dx1;
              
              % Any value in DXY2 which has a value of 2 is a common
              % derivative to both.
              dxyref = dxy2(:) == 2;
              % Can now get our mapping from DYI and DXI
              DerInds{Vcount,Icount}(dyi(dxyref),Fcount) = dxi(dxyref);
            end
          end
        end
      end
    end
  end
end
return
end

function DerInds = GetVertcatInds(Sizes,Vars)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;
NUMinput = length(Vars)-1;
NUMloops = size(Sizes,2);
y = Vars{1};
yMrow = y.func.size(1); yNcol = y.func.size(2);
if isinf(yNcol); yNcol = 1; end
yOverRef = zeros(yMrow,yNcol); yOverRef(:) = 1:yMrow*yNcol;
Inputs = Vars(2:end);

DerInds = cell(NUMvod,NUMinput);
[~,DPflag] = cadafuncname(y.id);

if DPflag
  yRowSizes = sum(Sizes(2:end,:),1);
  yColSizes = Sizes(1,:);
  iRowSizes = Sizes(2:end,:);
  iRowSums = zeros(NUMinput,NUMloops);
  for Icount = 1:NUMinput
    iRowSums(Icount,:) = sum(iRowSizes(1:Icount,:),1);
  end
  
  % See if Y is changing Col Size
  yColsTemp = nonzeros(yColSizes);
  if sum(yColsTemp ~= yColsTemp(1)); yColDep = 1; else yColDep = 0; end
  % If y is changing col size then everything is changing size
  
  if ~yColDep
    % if y cols arent changing, then check to see what is changing.
    yRowsTemp = nonzeros(yRowSizes);
    if sum(yRowsTemp ~= yRowsTemp(1)); yRowDep = 1; else yRowDep = 0; end
    
    iRowsDep = zeros(NUMinput,1);
    for Icount = 1:NUMinput
      iRowsTemp = nonzeros(Sizes(Icount+1,:));
      if sum(iRowsTemp ~= iRowsTemp(1)); iRowsDep(Icount) = 1; end
    end
  end
  
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      yrows = y.deriv(Vcount).nzlocs(:,1); ycols = y.deriv(Vcount).nzlocs(:,2);
      nzy = length(yrows);
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % Overmapped DY with linear index of DY as the values.
      for Icount = 1:NUMinput
        x = Inputs{Icount};
        % X = Input(I);
        if ~isempty(x.deriv(Vcount).nzlocs);
          DerInds{Vcount,Icount} = zeros(nzy,NUMloops);
          xMrow = x.func.size(1); xNcol = x.func.size(2);
          if isinf(xNcol); xNcol = 1; end
          xrows = x.deriv(Vcount).nzlocs(:,1); xcols = x.deriv(Vcount).nzlocs(:,2);
          nzx = length(xrows);
          dxoveri = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
          % Overmapped DX with linear index of DX as the values.
          xOverRef = zeros(xMrow,xNcol); xOverRef(:) = 1:xMrow*xNcol;
          for Fcount = 1:NUMloops
            if yRowSizes(Fcount)
              % This loop has a horzcat in it.
              if yColDep || (yRowDep && iRowsDep(Icount))
                % Both x and y are changing size
                yref = yOverRef(1:yRowSizes(Fcount),1:yColSizes(Fcount));
                % YREF - Function indices corresponding to Y on this loop.
                xref = xOverRef(1:iRowSizes(Icount,Fcount),1:yColSizes(Fcount));
                % XREF - Function indices corresponding to X on this loop.
              elseif yRowDep
                % Y is changing size, this input is not.
                xref = xOverRef;
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef(1:yRowSizes(Fcount),1:yColSizes(Fcount));
                % YREF - Function indices corresponding to Y on this loop.
              elseif iRowsDep(Icount)
                % Y is not changing size, but this input is.
                xref = xOverRef(1:iRowSizes(Icount,Fcount),1:yColSizes(Fcount));
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef;
                % YREF - Function indices corresponding to Y on this loop.
              else
                % Neither Y nor X is changing size
                xref = xOverRef;
                % XREF - Function indices corresponding to X on this loop.
                yref = yOverRef;
                % YREF - Function indices corresponding to Y on this loop.
              end
              dxi = dxoveri(xref(:),:);
              % DXI - Section of Overmapped DX corresponding to this loop -
              % with the linear index of DX as the values
              
              % Need to look at the part of DYI that DX is getting mapped
              % into
              if Icount == 1
                XintoYref = yref(1:iRowSizes(1,Fcount),:);
              else
                XintoYref = yref(iRowSums(Icount-1,Fcount)+1:iRowSums(Icount,Fcount),:);
              end
              % XintoYref - Function indices of Overmapped Y which the X in
              % this loop is getting mapped into.
              dyi = dyoveri(XintoYref(:),:);
              % DYI - Section of Overmapped DY which corresponds to where
              % DXI is getting mapped into in this loop.
              
              % Now have DXI and DYI, need to see what derivatives are
              % common to both of them.
              [yrows,ycols] = find(dyi); [xrows,xcols] = find(dxi);
              dy1 = sparse(yrows,ycols,ones(size(yrows)),size(dyi,1),size(dyi,2));
              dx1 = sparse(xrows,xcols,ones(size(xrows)),size(dxi,1),size(dxi,2));
              dxy2 = dy1+dx1;
              
              % Any value in DXY2 which has a value of 2 is a common
              % derivative to both.
              dxyref = dxy2(:) == 2;
              % Can now get our mapping from DYI and DXI
              DerInds{Vcount,Icount}(dyi(dxyref),Fcount) = dxi(dxyref);
            end
          end
        end
      end
    end
  end
end
return
end

function DerInds = GetTransposeInds(Sizes,y,x)
% Sizes are sizes of Y for each loop
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;

nzRef       = logical(sum(Sizes,1));
nnzLoops    = nnz(nzRef);
nLoops      = size(Sizes,2);
y.func.size(isinf(y.func.size)) = 1;
x.func.size(isinf(x.func.size)) = 1;
yMrow       = y.func.size(1); 
yNcol       = y.func.size(2);
xMrow       = x.func.size(1); 
xNcol       = x.func.size(2);
yOverRef    = zeros(yMrow,yNcol); 
yOverRef(:) = 1:yMrow*yNcol;
xOverRef    = zeros(xMrow,xNcol); 
xOverRef(:) = 1:xMrow*xNcol;

[~,DPflag] = cadafuncname(y.id);
DerInds = cell(NUMvod,1);
if DPflag
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      yrows = y.deriv(Vcount).nzlocs(:,1);
      ycols = y.deriv(Vcount).nzlocs(:,2); nzy = length(yrows);
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % OverMapped DY with linear indices as values
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2); nzx = length(xrows);
      dxoveri = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
      % OverMapped DX with linear indices as values
      DerIndsTemp = zeros(nzy,nnzLoops);
      for Fcount = 1:nnzLoops
        % Get the function locations of X for this loop
        yM = Sizes(1,Fcount); yN = Sizes(2,Fcount);
        xref = xOverRef(1:yN,1:yM);
        % Get the function locations of Y for this loop
        yref = yOverRef(1:yM,1:yN);
        
        % Reference Off what part of dyoveri is associated with this Y
        dyi = dyoveri(yref(:),:);
        % Reference off what part of dxoveri is associated with this Y
        dxi = dxoveri(xref(:),:);
        
        % Need to "transpose" the function indices of dxi
        yref = zeros(yN,yM);
        yref(:) = 1:yM*yN;
        yref = yref.';
        dxi = dxi(yref(:),:);
        
        % Find what derivatives are common to both DYI and DXI
        [yrows,ycols] = find(dyi); [xrows,xcols] = find(dxi);
        dy1 = sparse(yrows,ycols,ones(size(yrows)),size(dyi,1),size(dyi,2));
        dx1 = sparse(xrows,xcols,ones(size(xrows)),size(dxi,1),size(dxi,2));
        dxy2 = dy1+dx1;
        
        % Any Value in DXY2 which has a value of 2 is a common derivative to
        % both.
        dxyref = dxy2(:)==2;
        % Can now get our Linear Mapping from DXI and DYI as they have the
        % linear indices of the X and Y derivatives
        DerIndsTemp(dyi(dxyref),Fcount) = dxi(dxyref);
      end
      DerInds{Vcount}          = zeros(nzy,nLoops);
      DerInds{Vcount}(:,nzRef) = DerIndsTemp;
    end
  end
end
return
end

function DerInds = GetRepmatInds(Sizes,y,x)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;

nzRef    = logical(sum(Sizes,1));
nnzLoops = nnz(nzRef);
nLoops   = size(Sizes,2);
ySizes   = Sizes(1:2,nzRef);
xSizes   = Sizes(3:4,nzRef);
repM     = ySizes(1,:)./xSizes(1,:);
repN     = ySizes(2,:)./xSizes(2,:);
y.func.size(isinf(y.func.size)) = 1;
x.func.size(isinf(x.func.size)) = 1;
yMrow = y.func.size(1); 
yNcol = y.func.size(2);
xMrow = x.func.size(1); 
xNcol = x.func.size(2);
yOverRef    = zeros(yMrow,yNcol); 
yOverRef(:) = 1:yMrow*yNcol;
xOverRef    = zeros(xMrow,xNcol); 
xOverRef(:) = 1:xMrow*xNcol;
[~,DPflag] = cadafuncname(y.id);
DerInds = cell(NUMvod,1);
if DPflag
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      yrows = y.deriv(Vcount).nzlocs(:,1);
      ycols = y.deriv(Vcount).nzlocs(:,2); nzy = length(yrows);
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % OverMapped DY with linear indices as values
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2); nzx = length(xrows);
      dxoveri = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
      % OverMapped DX with linear indices as values
      
      DerIndsTemp = zeros(nzy,nnzLoops);
      for Fcount = 1:nnzLoops
        % Get the function locations of X for this loop
        xref = xOverRef(1:xSizes(1,Fcount),1:xSizes(2,Fcount));
        % Get the function locations of Y for this loop
        yref = yOverRef(1:ySizes(1,Fcount),1:ySizes(2,Fcount));
        % Get the Mapping of X when it is RepMatted.
        xRepRef = repmat(xref,repM(Fcount),repN(Fcount));
        
        % Reference Off what part of dyoveri is associated with this Y
        dyi = dyoveri(yref(:),:);
        % Reference off what part of dxoveri is associated with this Y
        dxi = dxoveri(xRepRef(:),:);
        
        % Find what derivatives are common to both DYI and DXI
        [yrows,ycols] = find(dyi); 
        [xrows,xcols] = find(dxi);
        dy1  = sparse(yrows,ycols,ones(size(yrows)),size(dyi,1),size(dyi,2));
        dx1  = sparse(xrows,xcols,ones(size(xrows)),size(dxi,1),size(dxi,2));
        dxy2 = dy1+dx1;
        
        % Any Value in DXY2 which has a value of 2 is a common derivative to
        % both.
        dxyref = dxy2(:)==2;
        % Can now get our Linear Mapping from DXI and DYI as they have the
        % linear indices of the X and Y derivatives
        DerIndsTemp(dyi(dxyref),Fcount) = dxi(dxyref);
      end
      DerInds{Vcount}          = zeros(nzy,nLoops);
      DerInds{Vcount}(:,nzRef) = DerIndsTemp;
    end
  end
end
return
end

function DerInds = GetReshapeInds(Sizes,y,x)
global ADIGATOR
NUMvod = ADIGATOR.NVAROFDIFF;

nzRef       = logical(sum(Sizes,1));
nnzLoops    = nnz(nzRef);
nLoops      = size(Sizes,2);
ySizes      = Sizes(1:2,nzRef);
xSizes      = Sizes(3:4,nzRef);
y.func.size(isinf(y.func.size)) = 1;
x.func.size(isinf(x.func.size)) = 1;
yMrow       = y.func.size(1); 
yNcol       = y.func.size(2);
xMrow       = x.func.size(1); 
xNcol       = x.func.size(2);
yOverRef    = zeros(yMrow,yNcol); 
yOverRef(:) = 1:yMrow*yNcol;
xOverRef    = zeros(xMrow,xNcol); 
xOverRef(:) = 1:xMrow*xNcol;

[~,DPflag] = cadafuncname(y.id);
DerInds = cell(NUMvod,1);
if DPflag
  for Vcount = 1:NUMvod
    if ~isempty(y.deriv(Vcount).nzlocs)
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      yrows = y.deriv(Vcount).nzlocs(:,1);
      ycols = y.deriv(Vcount).nzlocs(:,2); nzy = length(yrows);
      dyoveri = sparse(yrows,ycols,1:nzy,yMrow*yNcol,nv);
      % OverMapped DY with linear indices as values
      xrows   = x.deriv(Vcount).nzlocs(:,1);
      xcols   = x.deriv(Vcount).nzlocs(:,2); 
      nzx     = length(xrows);
      dxoveri = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
      % OverMapped DX with linear indices as values
      DerIndsTemp = zeros(nzy,nnzLoops);
      for Fcount = 1:nnzLoops
        % Get the function locations of X for this loop
        xref = xOverRef(1:xSizes(1,Fcount),1:xSizes(2,Fcount));
        % Get the function locations of Y for this loop
        yref = yOverRef(1:ySizes(1,Fcount),1:ySizes(2,Fcount));
        
        % Reference Off what part of dyoveri is associated with this Y
        dyi = dyoveri(yref(:),:);
        % Reference off what part of dxoveri is associated with this Y
        dxi = dxoveri(xref(:),:);
        
        % Find what derivatives are common to both DYI and DXI
        [yrows,ycols] = find(dyi); 
        [xrows,xcols] = find(dxi);
        dy1  = sparse(yrows,ycols,ones(size(yrows)),size(dyi,1),size(dyi,2));
        dx1  = sparse(xrows,xcols,ones(size(xrows)),size(dxi,1),size(dxi,2));
        dxy2 = dy1+dx1;
        
        % Any Value in DXY2 which has a value of 2 is a common derivative to
        % both.
        dxyref = dxy2(:)==2;
        % Can now get our Linear Mapping from DXI and DYI as they have the
        % linear indices of the X and Y derivatives
        DerIndsTemp(dyi(dxyref),Fcount) = dxi(dxyref);
      end
      DerInds{Vcount}          = zeros(nzy,nLoops);
      DerInds{Vcount}(:,nzRef) = DerIndsTemp;
    end
  end
end
return
end

function DataNameCount = AssignReferenceNames(Data,DATADEP,DATALOCS,...
  FORLENGTHS,FunStr,AsgnStr,AsgnLoc,AsgnFlag,DataNameCount)

% Used to Assign Reference Names for Indexing of Special Functions within a
% FOR loop.
% Inputs:
%   Data -  The data which is to be referenced for sizes or derivative
%           indexing.
%   DATADEP - Ones/Zeros of length of how many FOR loops this special
%             function is embedded within.
%   DATALOCS -  The Special Function Count and the For Count for each FOR 
%               loop this special function is embedded within.
%   FORLENGTHS -  Length of each FOR loop that this special function is
%                 embedded within.
%   FunStr  - Tells what special function this is (ex. '(FunStr)',
%             'SUBSASGN'
%   AsgnStr - Either 'INDICES' or 'SIZE' - tells where to assign these
%             indice names to.
%   AsgnLoc - Tells what row index of cells to place stuff into.
%   AsgnFlag - Tells whether or not to Assign to Cells 1:3 or 4:6
% We assign then as
% if ~AsgnFlag
%   ADIGATORFORDATA(DATALOCS(1,i)).(FunStr)(DATALOCS(2,i).(AsgnStr){AsgnLoc,1:3}
% else
%   ADIGATORFORDATA(DATALOCS(1,i)).(FunStr)(DATALOCS(2,i).(AsgnStr){AsgnLoc,4:6}
% Where 'FunStr' and 'AsgnStr' are structure field references and 'i' is
% is an index across all loops within which this particular function is
% embedded.

global ADIGATOR ADIGATORDATA ADIGATORFORDATA 
DepForLengths = FORLENGTHS(DATADEP);
DepDataLocs   = DATALOCS(1,DATADEP);

if DATADEP(1); ADIGATOR.FUNDEP = 1; end

if ~isempty(DepDataLocs)
  Data = reshape(Data,numel(Data)/DepForLengths(1),DepForLengths(1));
  CountName = ADIGATORFORDATA(DepDataLocs(1)).COUNTNAME;
  DepDataLocs(1) = []; DepForLengths(1) = [];
end

% ------------------------ Assign Data Name ----------------------------- %
if nnz(Data) < numel(Data); 
  SpFlag = 1; Data = sparse(Data);
else
  SpFlag = 0;
end
INDEXCOUNT = ADIGATORDATA.INDEXCOUNT + 1;
ADIGATORDATA.INDEXCOUNT = INDEXCOUNT;
ADIGATORDATA.INDICES.(sprintf('Index%1.0d',INDEXCOUNT)) = Data;
DataName = sprintf('Gator%1.0dIndices.Index%1.0d',...
  ADIGATOR.DERNUMBER,INDEXCOUNT);

% -------------- Middle FOR loops - reshaping stuff --------------------- %
for Fcount = 1:length(FORLENGTHS)-1
  if DATADEP(Fcount)
    % The Data is Dependent upon this FOR loop - Reshaping and Referencing
    % will take place before the next FOR loop runs.
    if ~isempty(DepForLengths)
      % There is another dependent loop past this one - This loop
      % references off of the Previous DataName, then reshapes this data
      % name.
      Data = reshape(Data(:,1),numel(Data(:,1))/DepForLengths(1),DepForLengths(1));
      DepForLengths(1) = [];
      % Make RHS that gets printed in this loop.
      ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
        {AsgnLoc,AsgnFlag+2} = ...
        sprintf(['reshape(',DataName,'(:,',CountName,'),%1.0d,%1.0d)'],...
        size(Data,1),size(Data,2));
      % Get Count Name of this FOR loop.
      CountName = ADIGATORFORDATA(DepDataLocs(1)).COUNTNAME;
      DepDataLocs(1) = [];
    else
      % There is not another dependent loop past this one - Just reference
      % off a column of the Previous Dataname.
      ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
        {AsgnLoc,AsgnFlag+2} = [DataName,'(:,',CountName,')'];
      % Make RHS string that gets printed at beginning of this loop.
    end
    DataNameCount = DataNameCount+1;
    DataName = sprintf('cada%1.0dforindex%1.0d',ADIGATOR.DERNUMBER,...
      DataNameCount);
    ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
      {AsgnLoc,AsgnFlag+1} = DataName;
    % Assign this loops Data Name.
  end
end

% ------Innermost FOR loop - what the Overloaded Operation will see------ %
ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
  {AsgnLoc,AsgnFlag+1} = DataName;
if DATADEP(end)
  ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
    {AsgnLoc,AsgnFlag+3} = [DATALOCS(1,end),SpFlag];
else
  ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
    {AsgnLoc,AsgnFlag+3} = [0,SpFlag];
end
return
end

function DataNameCount = AssignReferenceNamesCat(Data,DATADEP,DATALOCS,...
  FORLENGTHS,FunStr,AsgnStr,AsgnLoc,DataNameCount)

% Used to Assign Reference Names for Indexing of Special Functions within a
% FOR loop.
% Inputs:
%   Data -  The data which is to be referenced for sizes or derivative
%           indexing.
%   DATADEP - Ones/Zeros of length of how many FOR loops this special
%             function is embedded within.
%   DATALOCS -  The Special Function Count and the For Count for each FOR 
%               loop this special function is embedded within.
%   FORLENGTHS -  Length of each FOR loop that this special function is
%                 embedded within.
%   FunStr  - Tells what special function this is (ex. '(FunStr)',
%             'SUBSASGN'
%   AsgnStr - Either 'INDICES' or 'SIZE' - tells where to assign these
%             indice names to.
%   AsgnLoc - Tells what row index of cells to place stuff into.
% We assign then as
% if ~AsgnFlag
%   ADIGATORFORDATA(DATALOCS(1,i)).(FunStr)(DATALOCS(2,i).(AsgnStr){AsgnLoc(1),AsgnLoc(2),1:3}
% else
%   ADIGATORFORDATA(DATALOCS(1,i)).(FunStr)(DATALOCS(2,i).(AsgnStr){AsgnLoc(1),AsgnLoc(2),4:6}
% Where 'FunStr' and 'AsgnStr' are structure field references and 'i' is
% is an index across all loops within which this particular function is
% embedded.

global ADIGATOR ADIGATORDATA ADIGATORFORDATA 
DepForLengths = FORLENGTHS(DATADEP);
DepDataLocs   = DATALOCS(1,DATADEP);

if DATADEP(1); ADIGATOR.FUNDEP = 1; end

if ~isempty(DepDataLocs)
  Data = reshape(Data,numel(Data)/DepForLengths(1),DepForLengths(1));
  CountName = ADIGATORFORDATA(DepDataLocs(1)).COUNTNAME;
  DepDataLocs(1) = []; DepForLengths(1) = [];
end

% ------------------------ Assign Data Name ----------------------------- %
if nnz(Data) < numel(Data); 
  SpFlag = 1; Data = sparse(Data);
else
  SpFlag = 0;
end
INDEXCOUNT = ADIGATORDATA.INDEXCOUNT + 1;
ADIGATORDATA.INDEXCOUNT = INDEXCOUNT;
ADIGATORDATA.INDICES.(sprintf('Index%1.0d',INDEXCOUNT)) = Data;
DataName = sprintf('Gator%1.0dIndices.Index%1.0d',...
  ADIGATOR.DERNUMBER,INDEXCOUNT);

% -------------- Middle FOR loops - reshaping stuff --------------------- %
for Fcount = 1:length(FORLENGTHS)-1
  if DATADEP(Fcount)
    % The Data is Dependent upon this FOR loop - Reshaping and Referencing
    % will take place before the next FOR loop runs.
    if ~isempty(DepForLengths)
      % There is another dependent loop past this one - This loop
      % references off of the Previous DataName, then reshapes this data
      % name.
      Data = reshape(Data(:,1),numel(Data(:,1))/DepForLengths(1),DepForLengths(1));
      DepForLengths(1) = [];
      % Make RHS that gets printed in this loop.
      ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
        {AsgnLoc(1),AsgnLoc(2),2} = ...
        sprintf(['reshape(',DataName,'(:,',CountName,'),%1.0d,%1.0d)'],...
        size(Data,1),size(Data,2));
      % Get Count Name of this FOR loop.
      CountName = ADIGATORFORDATA(DepDataLocs(1)).COUNTNAME;
      DepDataLocs(1) = [];
    else
      % There is not another dependent loop past this one - Just reference
      % off a column of the Previous Dataname.
      ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
        {AsgnLoc(1),AsgnLoc(2),2} = [DataName,'(:,',CountName,')'];
      % Make RHS string that gets printed at beginning of this loop.
    end
    DataNameCount = DataNameCount+1;
    DataName = sprintf('cada%1.0dforindex%1.0d',ADIGATOR.DERNUMBER,...
      DataNameCount);
    ADIGATORFORDATA(DATALOCS(1,Fcount+1)).(FunStr)(DATALOCS(2,Fcount+1)).(AsgnStr)...
      {AsgnLoc(1),AsgnLoc(2),1} = DataName;
    % Assign this loops Data Name.
  end
end

% ------Innermost FOR loop - what the Overloaded Operation will see------ %
ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
  {AsgnLoc(1),AsgnLoc(2),1} = DataName;
if DATADEP(end)
  ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
    {AsgnLoc(1),AsgnLoc(2),3} = [DATALOCS(1,end),SpFlag];
else
  ADIGATORFORDATA(DATALOCS(1,end)).(FunStr)(DATALOCS(2,end)).(AsgnStr)...
    {AsgnLoc(1),AsgnLoc(2),3} = [0,SpFlag];
end
return
end