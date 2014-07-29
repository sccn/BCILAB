function x = cadaForSubsindex(x,FunString)
% function y = cadaForSubsindex(x,FunString)
% -----I THINK THIS FUNCTION IS NOW OBSOLETE!!!!!!!!! --------------------
%
% This function is called when an overloaded object is being used as a
% reference off of a cell or structure. In order to get here, the function
% subsindex is overloaded and sets a flag. The adigatorVarAnalzyer then sees
% this flag and calls this function. Note: this only gets called if being
% called from within a FOR loop, otherwise we are fine. Here we only allow
% this to take place if 1. no other operations are performed on this line
% (i.e. y = s(i).field; is okay, y = x(s(i).field) is not okay), 2. the
% result of y = s(i).field has no derivatives. This is basically set up as
% a way to reference problem auxillery data stored within structure/cell
% arrays when in a FOR loop.
% --------------------- Input Information ------------------------------- %
% x:          this is the result of the reference that took place
% FunString:  this is the users string that was evaluated
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
OTHERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.OTHER + 1;
ADIGATORFORDATA(INNERLOC).COUNT.OTHER = OTHERCOUNT;
if ~ADIGATOR.RUNFLAG
  x = cadaEmptyEval(x);
  return
elseif ADIGATOR.RUNFLAG == 1
  % Check to make sure that the resulting variable has no derivatives
  % - adigatorVarAnalyzer will do the overmapping for us.
  x.id   = ADIGATOR.VARINFO.COUNT;
  for Vcount = 1:ADIGATOR.NVAROFDIFF
    if ~isempty(x.deriv(Vcount).nzlocs)
      error(['If using an overloaded object as a reference off of a cell/structure ',...
        'within a FOR loop (or sub-function), the result must not contain derivatives.',...
        'This ability is only designed for referencing off of auxillery data stored ',...
        'within cells/structures']);
    end
  end
  if ADIGATOR.EMPTYFLAG
    x.func.size = [0 0];
  end
  % Here in .OTHER(OTHERCOUNT).DATA we store the size of x, if 0 in either
  % element, then the size changes.
  if isempty(ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA)
    ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA = x.func.size;
  else
    oldRow = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1);
    oldCol = ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(2);
    if oldRow ~= x.func.size(1)
      ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1) = 0;
    end
    if oldCol ~= x.func.size(2)
      ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(2) = 0;
    end
  end
else
  % ------------------------ Printing Run ------------------------------- %
  fid     = ADIGATOR.PRINT.FID;
  indent  = ADIGATOR.PRINT.INDENT;
  xOld    = x;
  
  % Get Overmapped Output
  xOverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.COUNT,1);
  x        = ADIGATORVARIABLESTORAGE.OVERMAP{xOverLoc};
  x.id = ADIGATOR.VARINFO.COUNT;
  [funcstr, ~] = cadafuncname();
  x.func.name = funcstr;
  
  % Get what the variables reference name is
  refID          = ADIGATOR.SUBSINDEXFLAG;

  oldRefStr      = ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(refID,1)};
  [newRefStr,~]  = cadafuncname(refID);
  
  oldStringMatches = cell(6,1);
  oldStringMatches{1} = ['\(',oldRefStr,'\).'];
  oldStringMatches{2} = ['\(',oldRefStr,','];
  oldStringMatches{3} = [',',oldRefStr,'\).'];
  oldStringMatches{4} = ['\{',oldRefStr,'\}'];
  oldStringMatches{5} = ['\{',oldRefStr,','];
  oldStringMatches{6} = [',',oldRefStr,'\}'];
  
  newStringMatches = cell(6,1);
  newStringMatches{1} = ['(',newRefStr,').'];
  newStringMatches{2} = ['(',newRefStr,','];
  newStringMatches{3} = [',',newRefStr,').'];
  newStringMatches{4} = ['{',newRefStr,'}'];
  newStringMatches{5} = ['{',newRefStr,','];
  newStringMatches{6} = [',',newRefStr,'}'];
  
  EqLoc = strfind(FunString,'=');
  FunString = strtrim(FunString(EqLoc+1:end));
  if strcmp(FunString(end),';')
    FunString = strtrim(FunString(1:end-1));
  end
  
  [start,match] = regexp(FunString,oldStringMatches,'start','split');
  for I = 1:6
    if ~isempty(start{I})
      match = match{I};
      % match{1} is the stuff to the left of what we searched, match{2} is
      % the stuff to the right of what we searched.
      if ~isempty(match{2})
        oldLoc = strfind(xOld.func.name,match{2});
        RHSstr = [match{1},newStringMatches{I},xOld.func.name(oldLoc:end)];
      else
        RHSstr = [match{1},newStringMatches{I},match{2}];
      end
      TF1 = sprintf('cada%1.0ftf1',ADIGATOR.NVAROFDIFF);
      if ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA
        % Neither Row Size nor Col Size Changes
        fprintf(fid,[indent,funcstr,' = ',RHSstr,';\n']);
      else
        fprintf(fid,[indent,TF1,' = ',RHSstr,';\n']);
        fprintf(fid,[indent,funcstr,' = zeros(%1.0f,%1.0f);\n'],x.func.size(1),x.func.size(2));
        if ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(1)
          % Column Size Changes
          fprintf(fid,[indent,funcstr,'(:,1:size(',TF1,',2)) = ',TF1,';\n']);
        elseif ADIGATORFORDATA(INNERLOC).OTHER(OTHERCOUNT).DATA(2)
          % Row Size Changes
          fprintf(fid,[indent,funcstr,'(1:size(',TF1,',1),:) = ',TF1,';\n']);
        else
          % Row and Col Change Sizes
          fprintf(fid,[indent,funcstr,'(1:size(',TF1,',1),1:size(',TF1,',2)) = ',TF1,';\n']);
        end
      end
      break
    end
  end
end
ADIGATOR.VARINFO.LASTOCC(x.id,1) = x.id;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
if ~isa(x,'cada')
  x = class(x,'cada');
end
end