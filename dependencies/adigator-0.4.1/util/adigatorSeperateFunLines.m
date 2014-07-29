function [FunStr,NUMFunStr] = adigatorSeperateFunLines(FunStrFull)
%function [FunStr,NUMFunStr] = adigatorSeperateFunLines(FunStrFull)
% This function seperates function lines when multiple commands are given
% on the same line
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

if ~strcmp(FunStrFull(1),'%')  
  % Find comments at end of line
  commentloc = findcomments(FunStrFull);
  if commentloc
    CommentStr = strtrim(FunStrFull(commentloc:end));
    FunStrFull = strtrim(FunStrFull(1:commentloc-1));
  end
  % ---------------- Check for Multiple Commands on Single Line --------- %
  FunStrLocs = strfind(FunStrFull,';');
  doubleSemiCheck = diff(FunStrLocs) == 1;
  if any(doubleSemiCheck)
    doubleSemi = 1;
    FunStrLocs = FunStrLocs([~doubleSemiCheck true]);
  else
    doubleSemi = 0;
  end
  % Check and make sure ';' isnt in vertcat
  if ~isempty(FunStrLocs)
    SquareLocs1 = strfind(FunStrFull,'[');
    if ~isempty(SquareLocs1)
      SquareLocs2 = strfind(FunStrFull,']');
      NUMSQL = length(SquareLocs1);
      for Scount = 1:NUMSQL
        % Need to find where bracket with Scount closes
        if Scount < NUMSQL
          for Ccount = Scount:NUMSQL-1
            if SquareLocs2(Ccount) < SquareLocs1(Ccount+1)
              CloseLoc = Ccount;
              break
            end
            if Scount == NUMSQL-1
              CloseLoc = NUMSQL;
            end
          end
        else
          CloseLoc = Scount;
        end
        for S2count = 1:length(FunStrLocs)
          if FunStrLocs(S2count) > SquareLocs1(Scount) &&...
              FunStrLocs(S2count) < SquareLocs2(CloseLoc)
            FunStrLocs(S2count) = 0;
          end
        end
      end
      FunStrLocs = nonzeros(FunStrLocs);
    end
    SquareLocs1 = strfind(FunStrFull,'{');
    if ~isempty(SquareLocs1)
      SquareLocs2 = strfind(FunStrFull,'}');
      NUMSQL = length(SquareLocs1);
      for Scount = 1:NUMSQL
        % Need to find where bracket with Scount closes
        if Scount < NUMSQL
          for Ccount = Scount:NUMSQL-1
            if SquareLocs2(Ccount) < SquareLocs1(Ccount+1)
              CloseLoc = Ccount;
              break
            end
            if Scount == NUMSQL-1
              CloseLoc = NUMSQL;
            end
          end
        else
          CloseLoc = NUMSQL;
        end
        for S2count = 1:length(FunStrLocs)
          if FunStrLocs(S2count) > SquareLocs1(Scount) &&...
              FunStrLocs(S2count) < SquareLocs2(CloseLoc)
            FunStrLocs(S2count) = 0;
          end
        end
      end
      FunStrLocs = nonzeros(FunStrLocs);
    end
  end
  % ----------------seperate multiple commands-----------------------
  if isempty(FunStrLocs) ||...
      (length(FunStrLocs) == 1 && FunStrLocs == length(FunStrFull))
    % only one line
    FunStr = cell(1,1);
    FunStr{1} = FunStrFull;
    NUMFunStr = 1;
  elseif FunStrLocs(end) ~= length(FunStrFull)
    % 2 or more lines
    % last line either comment or unsuppressed
    NUMFunStr = length(FunStrLocs) + 1;
    FunStr = cell(NUMFunStr,1);
    FunStr{1} = FunStrFull(1:FunStrLocs(1));
    for FScount = 1:NUMFunStr-2
      FunStr{FScount+1} = strtrim(FunStrFull(FunStrLocs(FScount)+1:FunStrLocs(FScount+1)));
    end
    FunStr{NUMFunStr} = strtrim(FunStrFull(FunStrLocs(end)+1:end));
  else
    % 2 or more lines, last line suppressed
    NUMFunStr = length(FunStrLocs);
    FunStr = cell(NUMFunStr,1);
    FunStr{1} = FunStrFull(1:FunStrLocs(1));
    for FScount = 2:NUMFunStr
      FunStr{FScount} = strtrim(FunStrFull(FunStrLocs(FScount-1)+1:FunStrLocs(FScount)));
    end
  end
  % Check for double ; at end of line
  if doubleSemi
    FunStr = regexprep(FunStr,';;$',';');
  end
  if commentloc
    FunStr{end+1,1} = CommentStr;
  end
else
  NUMFunStr = 1;
  FunStr = cell(1);
  FunStr{1} = FunStrFull;
end
end
function comment = findcomments(str)
% Look for a comment at end of string
comment = 0;
commentlocs = strfind(str,'%');
if ~isempty(commentlocs)
  % Make sure no %'s are used to build a string
  leftstrlocs = regexp(str,'[=''\(,]\s*''');
  rightstrlocs = regexp(str,'''\s*[''\),]');
  % if these lengths arent equal its possible that they occur after the
  % comment
  if length(leftstrlocs) < length(rightstrlocs)
    rightstrlocs = rightstrlocs(1:length(leftstrlocs));
  elseif length(leftstrlocs) > length(rightstrlocs)
    leftstrlocs = leftsrlocs(1:length(rightstrlocs));
  end
  if ~isempty(leftstrlocs)
    for C = commentlocs
      if ~any(C > leftstrlocs & C < rightstrlocs)
        comment = C;
      end
    end
  else
    comment = commentlocs(1);
  end
end
end