% readcogdata() - Read CogniScan data
%
% Usage:
%   >> [eegdata,D,N,fs,gain,events] = 
%	readcogdata(filename,channel,gaincorrect,duration,offset,show,m);
%
% Input:
% 'filename'    - [string] full path to data file
% 'channel'     - channels to read [all channels]
% 'gaincorrect' - 1 - apply gain; 0 - do no apply gain [1] 
% 'duration'    - number of samples to read [all samples]
% 'offset'      - sample number to begin reading [1]
% 'm'           - output of calibrate.m to discretize event channel
% 'show'        - 1 - display progress; 0 - do not display progress [1]
%
% Output:
% 'eegdata'     - [channels x samples] EEG data
% 'D'           - number of channels in filename
% 'N'           - number of samples in data file 
% 'fs'          - sampling rate
% 'gain'        - gain
% 'events'      - discretized eegdata (useful for event channel)
%    
% Example:
% [m,s]= calibrate(fileevent); % e.g. fileevent=jumptest 
% [D,N,fs,gain] = readheader(filein);
% [eventchannel,D,N,fs,gain,discreteeventchannel]=readcogdata(filein,64,0,N,1,m,1);
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson and Lucas Parra
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


function [eegdata,D,N,fs,gain,events]=readcogdata(filename,channel,gaincorrect,duration,offset,show,m);

[D,N,fs,gain] = readheader(filename); % Load header data

if nargin<2, channel=1:D; end
if channel==0, channel=1:D; end

if nargin<4, duration=N; end
if duration==0, duration=N; end
if nargin<5,
    offset=0; 
end
if nargin>=5, 
    if isempty(offset), offset=1; end; 
    offset=offset-1; % Read file including offset sample 
end 
if nargin<3, gaincorrect=1; end
if nargin<6, show=1; end

fid=fopen([filename '.bin'],'r','b');

eegdata=[]; events=[];
if show,
    disp(['Reading ' filename '.bin']);
end

for offseti=1:length(offset),
    fseek(fid,2*D*offset(offseti),'bof');
    
    %blklen=50000;
    blklen=2^14;
    
    slices=floor(duration/blklen);
    
    for i=1:slices,
        slice = fread(fid,[D blklen],'int16');
        if gaincorrect, slice=gain.*slice; end
        
        eegdata=[eegdata slice(channel,:)];
        if show,
            fprintf(['\roffset: %d/%d - slice: %d/%d'], ...
                offseti,length(offset),i.*blklen,duration); 
        end
        
        if nargin==7,
            events = [events eventnum(slice(end,:),m)];    
        end
        
    end
    
    blklenlast=duration-(slices.*blklen);
    slice = fread(fid,[D blklenlast],'int16');
    if gaincorrect, slice=slice.*gain; end
    eegdata=[eegdata slice(channel,:)];
    
    if nargin==7,
        events = [events eventnum(slice(end,:),m)];    
    end
    
    if show,
        fprintf(['\roffset: %d/%d - slice: %d/%d'],offseti,length(offset),slices.*blklen+blklenlast,duration);
    end
    
end

if show, fprintf('\n'); end

function [D,N,fs,gain] = readheader(filename)

% [D,N,fs,gain] = readheader(file)

fid = fopen([filename '.hdr'],'r');
line = fgetl(fid);
line = fgetl(fid); D = sscanf(line,'%d');
line = fgetl(fid); N = sscanf(line,'%d');
line = fgetl(fid); fs = 1/sscanf(line,'%f');
line = fgetl(fid); gain = sscanf(line,'%f');
fclose(fid);
