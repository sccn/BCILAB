function [nch, npoints, freq, d] = ddt(filename)
% ddt(filename) Read data from a .ddt file
%
% [nch, npoints, freq, d] = ddt(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   nch - number of channels
%   npoints - number of data points for each channel
%   freq - A/D frequency
%   d - [nch npoints] data array 

if(nargin ~= 1)
   disp('1 input argument is required')
   return
end
nch = 0;
npoints = 0;

if(isempty(filename))
   [fname, pathname] = uigetfile('*.ddt', 'Select a ddt file');
	filename = strcat(pathname, fname);
end
fid = fopen(filename, 'r');
if(fid == -1)
	disp('cannot open file');
   return
end
disp(strcat('file = ', filename));
version = fread(fid, 1, 'int32');
dataoffset = fread(fid, 1, 'int32');
freq = fread(fid, 1, 'double');
nch = fread(fid, 1, 'int32');

fseek(fid, 0, 1);
fsize = ftell(fid);
frewind(fid);
fseek(fid, dataoffset, 0);
npoints = (fsize - dataoffset)/(nch*2);
d = fread(fid, [nch npoints], 'int16');

fclose(fid);
