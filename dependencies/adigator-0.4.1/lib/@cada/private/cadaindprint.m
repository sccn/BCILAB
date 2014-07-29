function varargout = cadaindprint(x,varargin)
% function cadaindprint(x,varargin)
% Overloaded index printing scheme called from overloaded CADA functions.
% x is the index which is to be printed, this function stores the index,
% then either gives back the reference which will be used in the evaluation
% of the created file to get to this index, or assigns the given input name
% to the index reference.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORDATA

if numel(x) == 1 && nargin == 1
  varargout{1} = sprintf('%1.0f',x);
  return
end
ADIGATORDATA.INDEXCOUNT = ADIGATORDATA.INDEXCOUNT+1;
INDEXCOUNT = ADIGATORDATA.INDEXCOUNT;
Ind = sprintf('Index%1.0d',INDEXCOUNT);
INDEXNAME = sprintf('Gator%1.0dIndices.Index%1.0d',ADIGATOR.DERNUMBER,INDEXCOUNT);
% if ~islogical(x) && isequal(x,ones(size(x)))
%   x = ones(numel(x),1);
% end
ADIGATORDATA.INDICES.(Ind) = x;
if nargin == 2
  fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,varargin{1},...
    ' = ',INDEXNAME,';\n']);
else
  varargout{1} = INDEXNAME;
end