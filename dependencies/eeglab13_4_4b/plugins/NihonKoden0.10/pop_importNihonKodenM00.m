% pop_importNihonKodenM00() - EEGLAB plugin for importing Nihon Koden .m00
%                             text-format data.
%
% Usage:
%   >> 
%
% Author: Makoto Miyakoshi SCCN, INC, UCSD
%
% See also: eegplugin_importNihonKodenM00() importNihonKodenM00()

% History:
% 07/17/2013 ver 1.0 by Makoto. Created for my collaboration with Eishi Asano.

% Copyright (C) 2013 Makoto Miyakoshi, SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG com] = pop_importNihonKodenM00(filename)

if isempty(filename)
    [hdrfile path] = uigetfile2('*.m00', 'Select Nihon Koden .m00 file - pop_importNihonKodenM00()');
    if nnz(hdrfile) == 0 && nnz(path) == 0
        error('Operation terminated by user input.')
    end
    filename = [path hdrfile];
else
    [path,hdrfile01,hdrfile02] = fileparts(filename);
    hdrfile = [hdrfile01 hdrfile02];
end
[output,srate,scale] = readNihonKodenM00(filename);
EEG = pop_importdata('dataformat','array','data', output,'srate',srate,'setname',hdrfile(1:end-4));
com = ['[EEG, com] = pop_importNihonKodenM00(''' filename ''');'];
return