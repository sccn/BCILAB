function s = adigatorCellStrucAsgn(svacount,s,b,AsgnString,varargin)
% This function is called whenever a user does a cell or structure
% assignment using a variable as an index - this is a workaround to the
% fact that subsasgn will not be called in this case.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

numInd = nargin - 4;
if isstruct(s); 
  structflag = 1;
  % Need to find the field we are assigning to
  s = orderfields(s);
  equalloc = strfind(AsgnString,'=');
  fieldloc = strfind(AsgnString,').');
  fieldloc = fieldloc(end);
  fieldname = strtrim(AsgnString(fieldloc+2:equalloc-1));
  nameloc  = strfind(AsgnString,'(');
  endloc   = fieldloc;
else
  structflag = 0;
  nameloc  = strfind(AsgnString,'{');
  endloc   = strfind(AsgnString,'}');
  endloc   = endloc(end);
end

nameloc = nameloc(end);
varname  = strtrim(AsgnString(1:nameloc-1));
if ADIGATOR.OPTIONS.PREALLOCATE
  % This determines our maximum size and where all we assign to
  oldSize = size(s);
  sdummy  = zeros(oldSize);
  if numInd == 1
    if length(varargin{1}) > 1
      error('Currently can only do single assignments to cell/structures')
    end
    sdummy(varargin{1}) = 1;
  elseif numInd == 2
    if length(varargin{1}) > 1 || length(varargin{2}) > 1
      error('Currently can only do single assignments to cell/structures')
    end
    sdummy(varargin{1},varargin{2}) = 1; 
  else
    error('Currently can only do cell/structure assignment in two dimensions')
  end
  [I,J] = find(sdummy);
  if structflag
    s(I,J).(fieldname) = b;
  else
    % Cell if not a structure
    s{I,J} = b;
  end
  if isempty(I)
    inds = zeros(0,2);
  else
    inds = [I,J];
  end

  if isempty(ADIGATOR.STRUCASGN(svacount).dummy)
    ADIGATOR.STRUCASGN(svacount).dummy = getDummy(s,[]);
    ADIGATOR.STRUCASGN(svacount).inds = inds;
  else
    ADIGATOR.STRUCASGN(svacount).dummy = ...
      getDummy(s,ADIGATOR.STRUCASGN(svacount).dummy);
    matchflag = 0;
    for i = 1:size(ADIGATOR.STRUCASGN(svacount).inds,1)
      if ADIGATOR.STRUCASGN(svacount).inds(i,:) == inds
        matchflag = 1;
        break
      end
    end
    if ~matchflag
      ADIGATOR.STRUCASGN(svacount).inds(i+1,:) = [I,J];
    end
  end
else
  % Pre-Allocate s
  sdummy = ADIGATOR.STRUCASGN(svacount).dummy;
  newsize = size(sdummy);
  oldsize = size(s);
  if structflag
    newfields = fieldnames(sdummy);
    % Fill in missing fields
    for Fcount = 1:length(newfields)
      if ~isfield(s,newfields{Fcount})
        for Icount = 1:newsize(1)
          for Jcount = 1:newsize(2)
            s(Icount,Jcount).(newfields{Fcount}) = ...
              sdummy(Icount,Jcount).(newfields{Fcount});
          end
        end
      end
    end
    %
    for Icount = newsize(1):-1:oldsize(1)+1
      for Jcount = 1:newsize(2)
        for Fcount = 1:length(newfields)
          s(Icount,Jcount).(newfields{Fcount}) = ...
              sdummy(Icount,Jcount).(newfields{Fcount});
        end
      end
    end
    for Icount = 1:newsize(1)
      for Jcount = newsize(2):-1:oldsize(2)+1
        for Fcount = 1:length(newfields)
          s(Icount,Jcount).(newfields{Fcount}) = ...
              sdummy(Icount,Jcount).(newfields{Fcount});
        end
      end
    end
    % Check to make sure an empty object isnt assigned there
    for Fcount = 1:length(newfields)
      for Icount = 1:newsize(1)
        for Jcount = 1:newsize(2)
          if (isa(s(Icount,Jcount).(newfields{Fcount}),'cada') || isempty(s(Icount,Jcount).(newfields{Fcount}))) && ...
              ~isempty(sdummy(Icount,Jcount).(newfields{Fcount}))
            s(Icount,Jcount).(newfields{Fcount}) = ...
              sdummy(Icount,Jcount).(newfields{Fcount});
          end
        end
      end
    end
  else
    for Icount = newsize(1):-1:oldsize(1)+1
      for Jcount = newsize(2):-1:oldsize(2)+1
          s{Icount,Jcount} = sdummy{Icount,Jcount};
      end
    end
  end

  if ~ADIGATOR.RUNFLAG
    % Empty Run - Determine if any calculations were performed to compute indices
    PreOpCount = ADIGATOR.PREOPCOUNT;
    if ADIGATOR.VARINFO.COUNT > PreOpCount
      ADIGATOR.VARINFO.NAMELOCS(PreOpCount:ADIGATOR.VARINFO.COUNT-1,2) = ...
        1:ADIGATOR.VARINFO.COUNT-1-PreOpCount;
    end
    
    for icount = 1:numInd
      % Assign Last Occurrence of Indices
      ind = varargin{icount};
      if isa(ind,'cada')
        ADIGATOR.VARINFO.LASTOCC(ind.id,1) = ADIGATOR.VARINFO.COUNT;
      end
    end    
  end
  if ~ADIGATOR.RUNFLAG || ...
      (ADIGATOR.RUNFLAG == 2 && ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL)
    % Either in Empty run or printing run of a FOR loop - do all
    % assignments
    asgnInds = ADIGATOR.STRUCASGN(svacount).inds;
    for i = 1:size(asgnInds,1)
      I = asgnInds(i,1); J = asgnInds(i,2);
      if structflag
        s(I,J).(fieldname) = b;
      else
        s{I,J} = b;
      end
    end
  else
    % Only do the current assignment.
    if numInd == 1
      if structflag
        s(varargin{1}.func.value).(fieldname) = b;
      else
        s{varargin{1}.func.value} = b;
      end
    else
      I = varargin{1};
      if isa(I,'cada')
        I = I.func.value;
      end
      J = varargin{2};
      if isa(J,'cada')
        J = J.func.value;
      end
      if structflag
        s(I,J).(fieldname) = b;
      else
        s{I,J} = b;
      end
    end
  end
  if ADIGATOR.RUNFLAG == 2
    % Need to tell VarAnalyzer what to print
    if numInd == 1
      AsgnString = [AsgnString(1:nameloc),varargin{1}.func.name,...
        AsgnString(endloc:end)];
    else
      if isa(varargin{1},'cada')
        Istr = varargin{1}.func.name;
      else
        Istr = sprintf('%1.0f',varargin{1});
      end
      if isa(varargin{2},'cada')
        Jstr = varargin{2}.func.name;
      else
        Jstr = sprintf('%1.0f',varargin{2});
      end
      AsgnString = [AsgnString(1:nameloc),Istr,',',Jstr,...
        AsgnString(endloc:end)];
    end
  end
  ADIGATOR.PREOPCOUNT = ADIGATOR.VARINFO.COUNT;
  
  % Send this through adigatorVarAnalyzer
  if structflag
    s4VA = s;
    for Fcount = 1:length(newfields)
      if ~strcmp(newfields{Fcount},fieldname)
        s4VA = rmfield(s4VA,newfields{Fcount});
      end
    end
    s4VA = adigatorVarAnalyzer(AsgnString,s4VA,varname,0);
    for Icount = 1:size(s,1)
      for Jcount = 1:size(s,2)
        s(Icount,Jcount).(fieldname) = s4VA(Icount,Jcount).(fieldname);
      end
    end
  else
    s = adigatorVarAnalyzer(AsgnString,s,varname,0);
  end 
end

end

function dummy = getDummy(new,old)
if isnumeric(new)
  if ~isnumeric(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to unroll']);
  end
  dummy = [];
elseif isstruct(new)
  newfields = fieldnames(new);
  if iscell(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to unroll']);
  elseif isstruct(old)
    newfields = unique([newfields;fieldnames(old)]);
  end
  newsize = size(new);
  oldsize = size(old);
  if length(newsize) > 2
    error('We restrict cells/structures to two dimensions');
  end
  if any(oldsize > newsize)
    error(['Cannot remove elements/fields of a structure/cell if you are ',...
      'referencing/assigning within a loop which you wish to unroll']);
  end
  dummy = new;
  for Scount = 1:length(newfields)
    dfield = newfields{Scount};
    if isnumeric(old) || ~isfield(old,dfield)
      for Icount = 1:newsize(1)
        for Jcount = 1:newsize(2)
          dummy(Icount,Jcount).(dfield) = getDummy(new(Icount,Jcount).(dfield),[]);
        end
      end
    elseif isfield(new,dfield)
      for Icount = 1:newsize(1)
        for Jcount = 1:newsize(2)
          if Icount <= oldsize(1) && Jcount <= oldsize(2)
            dummy(Icount,Jcount).(dfield) = getDummy(new(Icount,Jcount).(dfield),...
              old(Icount,Jcount).(dfield));
          else
            dummy(Icount,Jcount).(dfield) = getDummy(new(Icount,Jcount).(dfield),[]);
          end
        end
      end
    else
      error(['Cannot remove elements/fields of a structure/cell if you are ',...
      'referencing/assigning within a loop which you wish to unroll']);
    end
  end
elseif iscell(new)
  if isstruct(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to unroll']);
  end
  dummy = new;
  newsize = size(new);
  oldsize = size(old);
  if length(newsize) > 2
    error('We restrict cells/structures to two dimensions');
  end
  if any(oldsize > newsize)
    error(['Cannot remove elements/fields of a structure/cell if you are ',...
      'referencing/assigning within a loop which you wish to unroll']);
  end
  for Icount = 1:newsize(1)
    for Jcount = 1:newsize(2)
      if Icount <= oldsize(1) && Jcount <= oldsize(2)
        dummy{Icount,Jcount} = getDummy(new{Icount,Jcount},old{Icount,Jcount});
      else
        dummy{Icount,Jcount} = getDummy(new{Icount,Jcount},[]);
      end
    end
  end
else
   % Dont Know what this is - probably a string
  dummy = new;
end

end