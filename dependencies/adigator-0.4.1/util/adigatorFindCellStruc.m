function x = adigatorFindCellStruc(x,xStr)
%function x = adigatorFindCellStruc(x,xStr)
% This routine recursively calls itself to find numeric and/or overloaded
% objects which are embedded within cells and/or structures.
% --------------------------- Input Information ------------------------- %
% x:    the cell or structure
% xStr: the user defined string name of the cell/structure
% -------------------------- Output Information ------------------------- %
% x:    the cell or structure with all of its old numeric fields now
%       overloaded
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if isstruct(x)
  % X is a structure - look through its fields and see what is what.
  x = orderfields(x);
  xFieldNames = fieldnames(x);
  if length(size(x)) > 2
    error('Structures arrays may only be in two dimensions.')
  end
  for Srowcount = 1:size(x,1)
    for Scolcount = 1:size(x,2)
      if numel(x) == 1
        xRefStr = xStr;
      else
        xRefStr = sprintf([xStr,'(%1.0d,%1.0d)'],Srowcount,Scolcount);
        if ADIGATOR.RUNFLAG == 0
          if size(x,1) == 1
            Kstr = sprintf('%1.0d',Scolcount);
            % Check for previous linear indexing
            ADIGATOR.VARINFO.NAMES = regexprep(ADIGATOR.VARINFO.NAMES,['\<',xStr,'\(',Kstr,'\)'],xRefStr);
          elseif size(x,2) == 1
            Kstr = sprintf('%1.0d',Srowcount);
            % Check for previous linear indexing
            ADIGATOR.VARINFO.NAMES = regexprep(ADIGATOR.VARINFO.NAMES,['\<',xStr,'\(',Kstr,'\)'],xRefStr);
          end
        end
      end
      for Fcount = 1:length(xFieldNames)
        yField = xFieldNames{Fcount};
        y = x(Srowcount,Scolcount).(yField);
        yStr = [xRefStr,'.',yField];
        if iscell(y)||isstruct(y)
          y = adigatorFindCellStruc(y,yStr);
        elseif isnumeric(y)||isa(y,'cada')
          y = FoundNumOrOver(y,yStr);
        else
%           error(['ADiGator only allows arrays, cells, and structures.',...
%             ' Including those objects embedded within cells/structures.']);
        end
        x(Srowcount,Scolcount).(yField) = y;
      end
    end
  end
elseif iscell(x)
  % X is a cell - look through its contents.
  if length(size(x)) > 2
    error('Cell arrays may only be in two dimensions.')
  end
  for Icount = 1:size(x,1)
    Istr = sprintf('%1.0d',Icount);
    for Jcount = 1:size(x,2)
      Jstr = sprintf('%1.0d',Jcount);
      y = x{Icount,Jcount};
      yStr = [xStr,'{',Istr,',',Jstr,'}'];
      if ADIGATOR.RUNFLAG == 0
        if size(x,1) == 1
          Kstr = sprintf('%1.0d',Jcount);
          % Check for previous linear indexing
          ADIGATOR.VARINFO.NAMES = regexprep(ADIGATOR.VARINFO.NAMES,['\<',xStr,'\{',Kstr,'\}'],yStr);
        elseif size(x,2) == 1
          Kstr = sprintf('%1.0d',Icount);
          % Check for previous linear indexing
          ADIGATOR.VARINFO.NAMES = regexprep(ADIGATOR.VARINFO.NAMES,['\<',xStr,'\{',Kstr,'\}'],yStr);
        end
      end
      if iscell(y)||isstruct(y)
        y = adigatorFindCellStruc(y,yStr);
      elseif isnumeric(y)||isa(y,'cada')
        y = FoundNumOrOver(y,yStr);
      else
%         error(['ADiGator only allows arrays, cells, and structures.',...
%           ' Including those objects embedded within cells/structures.']);
      end
      x{Icount,Jcount} = y;
    end
  end
end
return
end

function x = FoundNumOrOver(x,xStr)
global ADIGATOR
if ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  % Both cases need to get their own Operation Count and Variable Name
  % Assigned.
  ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1) = ...
    ADIGATOR.VARINFO.COUNT;
  % Set VarName
  adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,xStr,0);
  if isnumeric(x)
    % x is Numeric
    x = adigatorMakeNumeric(x);
  else
    % x is overloaded
    ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
    if ADIGATOR.VARINFO.NAMELOCS(x.id,3) == inf
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = inf;
    elseif ADIGATOR.VARINFO.NAMELOCS(x.id,3) == -inf
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = -inf;
    end
    if ADIGATOR.VARINFO.NAMELOCS(x.id,2) < 0
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) = ADIGATOR.VARINFO.NAMELOCS(x.id,2);
    end
    x.id = ADIGATOR.VARINFO.COUNT;
  end
  
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if isnumeric(x)
    % x is Numeric
    x = adigatorMakeNumeric(x);
    x = cadaOverMap(x);
  else
    % x is Overloaded
    OverOp = x.id;
    x.id   = ADIGATOR.VARINFO.COUNT;
    x      = cadaOverMap(x,OverOp);
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  if isnumeric(x)
    % x is Numeric
    x = adigatorMakeNumeric(x);
    x = cadaOverMap(x);
  else
    % x is overloaded
    xOp   = ADIGATOR.VARINFO.COUNT;
    oldID = x.id;
    x.id  = xOp;

    if isinf(ADIGATOR.VARINFO.NAMELOCS(oldID,2))
      ADIGATOR.VARINFO.NAMELOCS(xOp,2) = Inf;
      ADIGATOR.VARINFO.NAMELOCS(xOp,3) = ADIGATOR.VARINFO.NAMELOCS(oldID,3);
    end
    x           = cadaOverMap(x);
    funcname    = cadafuncname(xOp);
    x.func.name = funcname;
    for Vcount = 1:ADIGATOR.NVAROFDIFF
      if ~isempty(x.deriv(Vcount).name)
        x.deriv(Vcount).name = cadadername(funcname,Vcount);
      end
    end
  end
  ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1) = ...
    ADIGATOR.VARINFO.COUNT;
end
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
return
end