function [y,s] = adigatorCellStrucRef(svacount,s,RefString,varargin)
% This function is called whenever a user does a cell or structure
% reference using a variable as an index - this is a workaround to the
% fact that subsasgn will not be called in this case.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
numInd = nargin - 3;
if isstruct(s); 
  s = orderfields(s);
  structflag = 1;
  % Need to find the field we are referencing
  fieldloc = strfind(RefString,').');
  fieldloc = fieldloc(end);
  fieldname = strtrim(RefString(fieldloc+2:end-1));
  if strcmp(fieldname(end),';')
    fieldname = strtrim(fieldname(1:end-1));
  end
  nameloc  = strfind(RefString,'(');
  endloc   = fieldloc;
else
  structflag = 0;
  nameloc  = strfind(RefString,'{');
  endloc   = strfind(RefString,'}');
  endloc   = endloc(end);
end
nameloc  = nameloc(end);
equalloc = strfind(RefString,'=');
equalloc = equalloc(1);
varname = strtrim(RefString(1:equalloc-1));
strucname = strtrim(RefString(equalloc+1:nameloc-1));
if ADIGATOR.OPTIONS.PREALLOCATE
  % --------------------------------------------------------------------- %
  %                          NUMERIC PREALLOCATION RUN                    %
  % --------------------------------------------------------------------- %
  % Determine what references are made
  sdummy = zeros(size(s));
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
    y = s(I,J).(fieldname);
  else
    % Cell if not a structure
    y = s{I,J};
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
  if ~ADIGATOR.RUNFLAG
    % Empty Run - Determine if any calculations were performed to compute indices
    PreOpCount = ADIGATOR.PREOPCOUNT;
    if ADIGATOR.VARINFO.COUNT > PreOpCount
      ADIGATOR.VARINFO.NAMELOCS(PreOpCount:ADIGATOR.VARINFO.COUNT-1,2) = ...
        1:ADIGATOR.VARINFO.COUNT-PreOpCount;
    end
    ADIGATOR.PREOPCOUNT = ADIGATOR.VARINFO.COUNT;
    for icount = 1:numInd
      % Assign Last Occurrence of Indices
      ind = varargin{icount};
      if isa(ind,'cada')
        ADIGATOR.VARINFO.LASTOCC(ind.id,1) = ADIGATOR.VARINFO.COUNT;
      end
    end  
  elseif ADIGATOR.RUNFLAG == 2
    % Need to tell VarAnalyzer what to print
    if numInd == 1
      RefString = [RefString(1:nameloc),varargin{1}.func.name,...
        RefString(endloc:end)];
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
      RefString = [RefString(1:nameloc),Istr,',',Jstr,...
        RefString(endloc:end)];
    end
  end
  if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ...
      size(ADIGATOR.STRUCASGN(svacount).inds,1) > 1
    % ------------------------------------------------------------------- %
    %            ROLLING A LOOP - DIFFERENT REFERENCE INDICES             %
    % ------------------------------------------------------------------- %
    % If we are rolling a loop, we have an issue that we are assigning
    % different objects to a variable depending upon loop iteration. Rather
    % than trying to determine how things change, we are going to act as if
    % the objects which we are referencing are assigned to from within the
    % loop. This will make the overmapping scheme pre-allocate all of the
    % objects to an overmapped size prior to the loop.
    
    
    % Pre-Allocate s
    AllInds = ADIGATOR.STRUCASGN(svacount).inds;
    inds = AllInds(1,:); % If not overmapping run just use the first ref.
    if ADIGATOR.RUNFLAG == 1 
      % Overmapping run - determine which index to reference
      if numInd == 1
        I = varargin{1};
        if isa(I,'cada')
          I = I.func.value;
        end
        sdummy = zeros(size(s));
        sdummy(I) = 1;
        [I,J] = find(sdummy);
      else
        I = varargin{1};
        if isa(I,'cada')
          I = I.func.value;
        end
        J = varargin{2};
        if isa(J,'cada')
          J = J.func.value;
        end
      end
      inds = [I,J];
    end
    
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
      for Icount = newsize(1):-1:oldsize(1)+1
        for Jcount = newsize(2):-1:oldsize(2)+1
          for Fcount = 1:length(newfields)
            s(Icount,Jcount).(newfields{Fcount}) = ...
              sdummy(Icount,Jcount).(newfields{Fcount});
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
    % Assign s
    s = adigatorVarAnalyzer([],s,strucname,0);
    if structflag
      % y that we are going to output.
      y = s(inds(1),inds(2)).(fieldname);
    else
      y = s{inds(1),inds(2)};
    end 

    % Assign y
    if isa(y,'cada') && ADIGATOR.RUNFLAG == 2
      % We don't want to send this to the VarAnalyzer, rather just print it
      % out here.
      
      newRefStr = strtrim(RefString(equalloc+1:end));
      if strcmp(newRefStr(end),';')
        newRefStr = strtrim(newRefStr(1:end-1));
      end
      
      funcname = cadafuncname(ADIGATOR.VARINFO.COUNT);
      if isinf(ADIGATOR.VARINFO.NAMELOCS(y.id,3))
        fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,funcname,' = ',...
          newRefStr,';\n']);
      else
        fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,varname,' = ',...
          newRefStr,';\n']);
      end
      y.id = ADIGATOR.VARINFO.COUNT;
      y.func.name = funcname;
      yderiv = y.deriv;
      for Vcount = 1:ADIGATOR.NVAROFDIFF
        if ~isempty(yderiv(Vcount).nzlocs)
          y.deriv(Vcount).name = cadadername(funcname,Vcount);
        end
      end
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    else
      y = adigatorVarAnalyzer(RefString,y,varname,0);
    end

    ADIGATOR.PREOPCOUNT = ADIGATOR.VARINFO.COUNT;
  else
    % ------------------------------------------------------------------- %
    %                       SINGLE REFERENCE INDEX                        %
    % ------------------------------------------------------------------- %
    % We can just do the reference for these and let adigatorVarAnalyzer
    % handle the rest.
    if ~ADIGATOR.RUNFLAG
      ind = ADIGATOR.STRUCASGN(svacount).inds(1,:);
      if structflag
        y = s(ind(1),ind(2)).(fieldname);
      else
        y = s{ind(1),ind(2)};
      end
    else
      if numInd == 1
        if isa(varargin{1},'cada')
          ind = varargin{1}.func.value;
        else
          ind = varargin{1};
        end
        if structflag
          y = s(ind).(fieldname);
        else
          y = s{ind};
        end
      else
        if isa(varargin{1},'cada')
          ind1 = varargin{1}.func.value;
        else
          ind1 = varargin{1};
        end
        if isa(varargin{2},'cada')
          ind2 = varargin{2}.func.value;
        else
          ind2 = varargin{2};
        end
        if structflag
          y = s(ind1,ind2).(fieldname);
        else
          y = s{ind1,ind2};
        end
      end
    end
    y = adigatorVarAnalyzer(RefString,y,varname,0);
  end
end

end

function dummy = getDummy(new,old)
if isnumeric(new)
  if ~isnumeric(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to keep rolled']);
  end
  dummy = [];
elseif isstruct(new)
  newfields = fieldnames(new);
  if iscell(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to keep rolled']);
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
      'referencing/assigning within a loop which you wish to keep rolled']);
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
      'referencing/assigning within a loop which you wish to keep rolled']);
    end
  end
elseif iscell(new)
  if isstruct(old)
    error(['If an element of a cell/structure is a cell/structure, ',...
      'and you are referencing/assigning to the parent cell/structure in a loop then ',...
      'the child element must stay a cell/structure in order to keep rolled']);
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