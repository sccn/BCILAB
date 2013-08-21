function timeinds = datetime2ind(datelist, timelist)
%DATETIME2IND      Converts date/time lists to serial time indices.
%   TIMEINDS = DATETIME2IND(DATELIST, TIMELIST) takes two N x 1 cell
%   arrays containing the dates and times (as strings) to be converted
%   and returns an N x 1 vector of (double) time indices as described by
%   'datenum'.
%
%   The DATELIST entries must be of the format 'mmddyy' and the TIMELIST
%   entries 'hhmm'.

% Argument checking.
if ((nargin < 2) || (length(datelist) ~= length(timelist)) || ...
	(~iscell(datelist)) || (~iscell(timelist)))
	error('Two cell arrays of equal length are required.');
end

if (isempty(datelist)),  timeinds = [];  return;   end

% Make copies of the dates and start times in a more user-friendly numeric matrix.
dateN = str2num(cat(1, datelist{:}));
timeN = str2num(cat(1, timelist{:}));

% Next, break dates and times into component pieces.
yearnums  = rem(dateN, 100);
monthnums = floor(dateN/10000);
daynums   = rem(floor(dateN/100), 100);
hournums  = floor(timeN/100);
minnums   = rem(timeN, 100);
secnums   = zeros(size(dateN,1), 1);

% Fix the yearnums (pivot on 2000)
latter = (yearnums <= 50);
yearnums(latter) = yearnums(latter) + 2000;
yearnums(~latter) = yearnums(~latter) + 1900;

% Use these pieces to compute the unique, sortable serial number for each record
timeinds = datenum(yearnums, monthnums, daynums, hournums, minnums, secnums);