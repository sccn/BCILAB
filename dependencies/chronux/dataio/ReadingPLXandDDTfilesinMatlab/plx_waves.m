function [n, npw, ts, wave] = plx_waves(filename, ch, u)
% plx_waves(filename, channel, unit): Read waveform data from a .plx file
%
% [n, npw, ts, wave] = plx_waves(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number
%   unit  - unit number (0- invalid, 1-4 valid)
% OUTPUT:
%   n - number of waveforms
%   npw - number of points in each waveform
%   ts - array of timestamps (in seconds) 
%   wave - array of waveforms [npw, n], raw a/d values

if(nargin ~= 3)
   disp('3 input arguments are required')
   return
end

n = 0;
npw = 0;
ts = 0;
wave = 0;

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

record = 0;
wave = zeros(npw, 1);
wf = zeros(npw, 1);

% read data records
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
      wf = fread(fid, [toread, 1], 'int16');
   end
   if toread > 0
   	if type == 1
         if channel == ch 
            if unit == u
 	       n = n + 1;
               ts(n) = timestamp/freq;
               wave(:, n) = wf;
            end
      	 end
   	end
   end
   
   record = record + 1;
   if feof(fid) == 1
      break
   end
   
end
disp(strcat('number of waveforms = ', num2str(n)));

fclose(fid);
