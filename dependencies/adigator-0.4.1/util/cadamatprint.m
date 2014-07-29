function varargout = cadamatprint(x,varargin)
% function cadamatprint(x,varargin)
% Overloaded matrix data printing scheme
% 
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORDATA

if numel(x) == 1 && nargin == 1
  varargout{1} = num2str(x,16);
  return
end
ADIGATORDATA.DATACOUNT = ADIGATORDATA.DATACOUNT+1;
DATACOUNT = ADIGATORDATA.DATACOUNT;
Ind = sprintf('Data%1.0d',DATACOUNT);
DATANAME = sprintf('Gator%1.0dData.Data%1.0d',ADIGATOR.DERNUMBER,DATACOUNT);

ADIGATORDATA.DATA.(Ind) = x;
if nargin == 2
  fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,varargin{1},...
    ' = ',DATANAME,';\n']);
else
  varargout{1} = DATANAME;
end