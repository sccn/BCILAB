function [n, ts, sv] = plx_event_ts(filename, ch)
% plx_event_ts(filename, channel) Read event timestamps from a .plx file
%
% [n, ts, sv] = plx_event_ts(filename, channel)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based external channel number
%             strobed channel has channel number 257  
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps 
%   sv - array of strobed event values (filled only if channel is 257)

if(nargin ~= 2)
   disp('2 input arguments are required')
   return
end

n = 0;
ts = 0;
sv = 0;

if(isempty(filename))
   [fname, pathname] = uigetfile('*.plx', 'Select a plx file');
	filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
	disp('cannot open file');
   return
end

disp(strcat('file = ', filename));

% read file header
header = fread(fid, 64, 'int32');
freq = header(35);  % frequency
ndsp = header(36);  % number of dsp channels
nevents = header(37); % number of external events
nslow = header(38);  % number of slow channels
npw = header(39);  % number of points in wave
npr = header(40);  % number of points before threshold
tscounts = fread(fid, [5, 130], 'int32');
wfcounts = fread(fid, [5, 130], 'int32');
evcounts = fread(fid, [1, 512], 'int32');

% skip variable headers 
fseek(fid, 1020*ndsp + 296*nevents + 296*nslow, 'cof');

% read the data
record = 0;
while feof(fid) == 0
	type = fread(fid, 1, 'int16');
	upperbyte = fread(fid, 1, 'int16');
	timestamp = fread(fid, 1, 'int32');
	channel = fread(fid, 1, 'int16');
	unit = fread(fid, 1, 'int16');
	nwf = fread(fid, 1, 'int16');
	nwords = fread(fid, 1, 'int16');
	toread = nwords;
	if toread > 0
	  wf = fread(fid, toread, 'int16');
	end
   	if type == 4
         if channel == ch 
 	        	n = n + 1;
         		ts(n) = timestamp;
				sv(n) = unit;
      	 end
   	end
   
   record = record + 1;
   if feof(fid) == 1
      break
   end
   
end
disp(strcat('number of timestamps = ', num2str(n)));

fclose(fid);
