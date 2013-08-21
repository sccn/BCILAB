function [tscounts, wfcounts, evcounts] = plx_info(filename, fullread)
% plx_info(filename, fullread) -- read and display .plx file info
%
% [tscounts, wfcounts] = plx_info(filename, fullread)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   fullread - if 0, reads only the file header
%              if 1, reads all the file
% OUTPUT:
%   tscounts - 5x130 array of timestamp counts
%      tscounts(i, j) is the number of timestamps for channel i, unit j
%   wfcounts - 5x130 array of waveform counts
%     wfcounts(i, j) is the number of waveforms for channel i, unit j
%   evcounts - 1x512 array of external event counts
%     evcounts(i) is the number of events for channel i

if(nargin ~= 2)
   disp('2 input arguments are required')
   return
end

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
header = fread(fid, 64, 'int32');
version = header(2);
freq = header(35);  % frequency
ndsp = header(36);  % number of dsp channels
nevents = header(37); % number of external events
nslow = header(38);  % number of slow channels
npw = header(39);  % number of points in wave
npr = header(40);  % number of points before threshold
disp(strcat('version = ', num2str(version)));
disp(strcat('frequency = ', num2str(freq)));
disp(strcat('number of DSP headers = ', num2str(ndsp)));
disp(strcat('number of Event headers = ', num2str(nevents)));
disp(strcat('number of A/D headers = ', num2str(nslow)));
tscounts = fread(fid, [5, 130], 'int32');
wfcounts = fread(fid, [5, 130], 'int32');
evcounts = fread(fid, [1, 512], 'int32');
if fullread > 0
   % reset counters
   tscounts = zeros(5, 130);
   wfcounts = zeros(5, 130);
   evcounts = zeros(1, 512);
   % skip variable headers
   fseek(fid, 1020*ndsp + 296*nevents + 296*nslow, 'cof');
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
      if type == 1
         tscounts(unit+1, channel+1) = tscounts(unit+1, channel+1) + 1;
         if toread > 0
            wfcounts(unit+1, channel+1) = wfcounts(unit+1, channel+1) + 1;
         end
      end
      if type == 4
         evcounts(channel+1) = evcounts(channel+1) + 1;
      end
      
   	record = record + 1;
	   if feof(fid) == 1
   	   break
	   end
	end
   disp(strcat('number of records = ', num2str(record)));
end
disp(' ');
disp(' Timestamps:');
disp(' ch unit  count');
for i=1:130
   for j=1:5
      if tscounts(j, i) > 0
         disp(sprintf('%3d %4d %6d', i-1, j-1, tscounts(j, i)));
      end
   end
end

disp(' ');
disp(' Waveforms:');
disp(' ch unit  count');
for i=1:130
   for j=1:5
      if wfcounts(j, i) > 0
         disp(sprintf('%3d %4d %6d', i-1, j-1, wfcounts(j, i)));
      end
   end
end

disp(' ');
disp(' Events:');
disp(' ch count');
for i=1:300
  if evcounts(i) > 0
     disp(sprintf('%3d %6d', i-1, evcounts(i)));
   end
end

disp(' ');
disp(' A/D channels:');
disp(' ch count');
for i=301:364
  if evcounts(i) > 0
     disp(sprintf('%3d %6d', i-301, evcounts(i)));
   end
end


fclose(fid);
