% readNihonKodenM00 - reads Nihon Koden .m00 text files.
%
% Usage  : [output,srate,scale] = readNihonKodenM00(input);
%
% Inputs : path to the Nihon Koden .m00 text file
% 
% Outputs: Matlab variable in channels x timepoints, sampling rate, and
%          amplitude scale. If 'scale' is 1, no need to rescale the data. If
%          'scale'~=1, rescale the data by the factor of 'scale'. 

% History:
% 07/17/2013 ver 1.0 by Makoto. Created for my collaboration with Eishi Asano.

% Copyright (C) 2013 Makoto Miyakoshi SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General PublIC LICense as published by
% the Free Software Foundation; either version 2 of the LICense, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General PublIC LICense for more details.
%
% You should have received a copy of the GNU General PublIC LICense
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function [output,srate,scale] = readNihonKodenM00(input)

% open the file
FID = fopen(input);

% skip the first two lines
firstLine  = fgetl(FID);
secondLine = fgetl(FID);

% parse header info in the first line
C = strsplit(firstLine, ' ');
timePnts = str2num(C{1,1}(12:end));             % 12 letters are 'TimePoints='
numChans = str2num(C{2,1}(10:end));             % 10 letters are 'Channels='
srate    = round(1000/str2num(C{4,1}(22:end))); % 22 letters are 'SamplingInterval[ms]='
scale    = str2num(C{5,1}(9:end));              % 9  letters are 'Bins/uV='

% parse data structure
tmpData = fscanf(FID, '%f');
fclose(FID);

% check consistency between header info and actual data size
if size(tmpData,1) ~= numChans*timePnts;
    warning('Data size does not match header information. Check consistency.')
end

% reshape the data for finalizing
output = reshape(tmpData, [numChans timePnts]);