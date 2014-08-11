function ind = subsindex(x)
% This will be called when someone uses a referece like structi = struct(i)
% Can't find these instances from the source code level.

% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.SUBSINDEXFLAG && ADIGATOR.FORINFO.FLAG
  error(['Cannot do multiple subsindexing with overloaded objects on the ',...
    'same line when in a FOR loop (or sub-function).']);
end
%if ADIGATOR.FORINFO.FLAG
%   if ~ADIGATOR.RUNFLAG
%     ind = 0;
%     ADIGATOR.SUBSINDEXFLAG = x.id;
%   elseif ADIGATOR.RUNFLAG == 1
%     if ADIGATOR.EMPTYFLAG
%       ind = 0;
%       ADIGATOR.SUBSINDEXFLAG = x.id;
%     else
%       if ~isempty(x.func.value)
%         ind = x.func.value-1;
%         ADIGATOR.SUBSINDEXFLAG = ind+1;
%       else
%         error('Cannot subsindex a strictly symbolic object')
%       end
%       if length(ind) > 1
%         error(['Can only single reference off of structure/cell array ',...
%           'with an overloaded object while in a FOR loop (or sub-function).']);
%       end
%     end
%   else
%     ind = 0;
%     ADIGATOR.SUBSINDEXFLAG = x.id;
%   end
if ADIGATOR.EMPTYFLAG  
  ind = 0;
elseif ~isempty(x.func.value)
  ind = x.func.value-1;
elseif ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ADIGATOR.PRINT.FLAG
  errStr = ['Attempting to reference a structure as ''mystructi = mystruct(i)'' OR',...
    ' assigning to a structure as ''mystruct(i) = mystructi'' ',...
    'cannot do either of these in a loop when reference index ''i'' changes - for ',...
    'references please use ''fieldi = mystruct(i).field'' or use cells as ',...
    '''mycelli = mycell{i}'' - for assignments please use ',...
    '''mystruct(i).field = fieldi'' or use cells as ''mycell{i} = mycelli'' ',...
    'only other alternative is to unroll loops'];
  error(errStr)
else
  error('Cannot subsindex a strictly symbolic object')
end