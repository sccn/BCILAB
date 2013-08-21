function [matchstrs,inds] = strmatchre(pattern, strs)
%STRMATCHRE        Select strings that match a regular expression.
%   MATCHSTRS = STRMATCHRE(PATTERN, STRS) returns the rows of the character
%   array or elements of the cell array of strings STRS that match the
%   regular expression PATTERN (see REGEXP for the syntax of PATTERN).
%   Matched strings must match the entire PATTERN -- so 'testing' matches
%   '.*ing' but not 'ing'.
%
%   If PATTERN is a cell array of regular expressions, MATCHSTRS contains
%   the strings from STRS that match at least one of the patterns.
%
%   [MATCHSTRS, INDS] = STRMATCHRE(PATTERN, STRS) also returns the indices
%   into STRS such that MATCHSTRS equals STRS(INDS,:) if STRS is a
%   character array or STRS{INDS} if STRS is a cell array.
%
%   Example:
%      D = dir;  filenames = {D.name};
%      STRMATCHRE(filenames, '.*\.m$')
%         returns the names of files in the current directory with 
%         extension .m
% 
%   See also STRMATCH, REGEXP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coerce strs to a 1 column cell array
origcell = iscell(strs);
if (~origcell),  strs = cellstr(strs);  end;
strs = strs(:);

% We'll need a list of string lengths later to make sure we matched all
% the way to the end
if (iscell(pattern)),  numpat = length(pattern);
else	               numpat = 1;   pattern = {pattern};
end
strlens = repmat(cellfun('length', strs), 1, numpat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Match Pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%
start = {};  finis = {};
for p = 1:numpat
	[sta,fin] = regexp(strs, pattern{p});
	start = cat(2,start,sta);  finis = cat(2,finis,fin);
end

%%%%%%%%%%%%%%%%%% Convert REGEXP output to indices %%%%%%%%%%%%%%%%%%
% Change failed matches from [] to 0 so we can convert to a matrix.
failed = (cellfun('length', start) ~= 1);  % (>1 match means failed too)
[start{failed}] = deal(0);   [finis{failed}] = deal(0);
start = cell2mat(start);     finis = cell2mat(finis);

% Only take strings that matched at least one pattern beginning to end
start(start ~= 1) = 0;       finis(finis ~= strlens) = 0;
inds = find(any(start,2) & any(finis,2));

%%%%%%%%%%%%%%%%%%%%%%%%% Construct outputs %%%%%%%%%%%%%%%%%%%%%%%%%%
matchstrs = strs(inds);
if (~origcell), matchstrs = char(matchstrs); end;	
if (nargout == 1), clear inds; end;