function [adfreq, n, ts, fn, ad] = plx_ad(filename, ch)
% plx_ad(filename, channel): Read a/d data from a .plx file
%
% [adfreq, n, ts, fn, ad] = plx_ad(filename, ch)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 0 - based channel number
%
%           a/d data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values stored in the vector ad. 
% OUTPUT:
%   n - total number of data points 
%   ts - array of fragment timestamps (one timestamp for fragment, in seconds)
%   fn - number of data points in each fragment
%   ad - array of raw a/d values

if(nargin ~= 2)
   disp('2 input arguments are required')
   return
end

i = 0;
n = 0;
ts = 0;
fn = 0;
ad = 0;

if(isempty(filename))
   [fname, pathname] = uigetfile('*.plx', 'Select a plx file');
	filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
	disp('cannot open file');
   return
end

% calculate file size
fseek(fid, 0, 'eof');
fsize = ftell(fid);
fseek(fid, 0, 'bof');


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

% A/D counts are stored in evcounts (301, 302, etc.)
count = 0;
if evcounts(301+ch) > 0
	count = evcounts(301+ch);
	ad = 1:count;
end

% skip DSP and Event headers
fseek(fid, 1020*ndsp + 296*nevents, 'cof');

% read one A/D header and get the frequency
adheader = fread(fid, 74, 'int32');
adfreq = adheader(10);

% skip all other a/d headers
fseek(fid, 296*(nslow-1), 'cof');

record = 0;

wf = zeros(1, npw);
adpos = 1;

while feof(fid) == 0
   type = fread(fid, 1, 'int16');
	upperbyte = fread(fid, 1, 'int16');
	timestamp = fread(fid, 1, 'int32');
	channel = fread(fid, 1, 'int16');
   unit = fread(fid, 1, 'int16');
   nwf = fread(fid, 1, 'int16');
   nwords = fread(fid, 1, 'int16');
   if nwords > 0
      wf = fread(fid, [1 nwords], 'int16');
   end
   if nwords > 0
      if type == 5
         if channel == ch 
        	i = i + 1;
			n = n + nwords;
         	ts(i) = timestamp/freq;
			fn(i) = nwords;
			if count > 0
			    if adpos+nwords-1 <= count
				    ad(adpos:adpos+nwords-1) = wf(1:nwords);
					adpos = adpos + nwords;
			    else
					for i=1:nwords
						ad(adpos) = wf(i); adpos = adpos + 1;
					end
			    end
			else
				ad = [ad wf(1, 1:nwords)];
			end
      	 end
      end
   end
   
   record = record + 1;
   if mod(record, 1000) == 0
       disp(sprintf('records %d points %d (%.1f%%)', record, n, 100*ftell(fid)/fsize));
   end

   if feof(fid) == 1
      break
   end
   
end

if adpos-1 < count
   ad = ad(1:adpos-1);
end

disp(strcat('number of data points = ', num2str(n)));

fclose(fid);
