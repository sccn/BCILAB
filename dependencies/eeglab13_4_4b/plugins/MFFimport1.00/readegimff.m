% readegimff() - read header, event info and data from EGI Meta File Format
% (MFF) recording.
%
% Usage:
%   >> [head,evt,data] = readegimff(filename,dtype,firstsample,lastsample)
%
% Inputs:
%   filename - string with full path to file
%   dtype    - type of data to read (EEG or PIB)
%   firstsample - index of first sample to read from data
%   lastsample - index of last sample to read from data
%
% Outputs:
%   head - data header structure containing the following fields
%     samp_rate: sampling rate
%         nchan: number of channels
%       samples: number of samples
%      segments: number of data segments
%      segsamps: number of samples/segment
%
%   evt - 1xN event structure array, with each element containing the
%   following fields
%          type: type of event
%        sample: sample point of event
%         value: value of event
%        offset: offset of event
%      duration: duration of event
%     timestamp: timestamp of event
%         codes: MX2 cell array of key code-value pairs
%
%   data - nchan X nsamples (continuous data) or nchan x nsamples x nsegments
%   (epoched data) double array
%
% See also:
%   eeglab(), pop_readegimff(), pop_readegi()
%
% Copyright (C) 2011 Srivas Chennu, University of Cambridge, srivas@gmail.com
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


function [head,evt,data] = readegimff(filename,dtype,firstsample,lastsample)

%%%% size of data chunk to read at one go. %%%%
CHUNKSIZE = 500000;

head = [];
evt = [];
data = [];

if nargin < 1
    help readegimffhdr;
    return;
end

if nargin < 2
    dtype = 0;
else
    if strcmp(dtype,'PIB')
        dtype = 1;
    else
        dtype = 0;
    end
end
    
mffjarpath = which('MFF-1.0.d0004.jar');
if isempty(mffjarpath)
    error('MFF-1.0.d0004.jar not found in path! Make sure the mffreader directory is in your MATLAB path');
end

javaaddpath(mffjarpath);


fprintf('Reading header.\n');
try
    mffhdr = read_mff_header(filename,dtype);
catch err
    if strcmp(err.identifier,'MATLAB:Java:GenericException') && ...
            ~isempty(strfind(err.message,'java.lang.OutOfMemoryError'))
        error(['%s\n\nreadegimff: Increase the size of your MATLAB Java heap space!\n' ...
            'Find out how at http://www.mathworks.co.uk/support/solutions/en/data/1-18I2C/'],...
            err.message);
    else
        rethrow(err);
    end
end

fprintf('Reading events.\n');
try
    evt = read_mff_event(filename);
catch err
    if strcmp(err.identifier,'MATLAB:Java:GenericException') && ...
            ~isempty(strfind(err.message,'java.lang.OutOfMemoryError'))
        error(['%s\n\nreadegimff: Increase the size of your MATLAB Java heap space!\n' ...
            'Find out how at http://www.mathworks.co.uk/support/solutions/en/data/1-18I2C/'],...
            err.message);
    else
        rethrow(err);
    end
end

fprintf('Reading data.\n');

if mffhdr.nTrials == 1
    if isempty(firstsample)
        firstsample = 1;
    elseif firstsample < 1
        error('%s: first sample to read must be >= 1.',mfilename);
    end
    
    if isempty(lastsample)
        lastsample = mffhdr.nSamples;
    elseif lastsample > mffhdr.nSamples
        error('%s: last sample to read must be <= %d.',mfilename,mffhdr.nSamples);
    end
    
    nSampRead = lastsample-firstsample+1;
    data = zeros(mffhdr.nChans,nSampRead);

    chunkStart = 1;
    while firstsample+chunkStart <= lastsample
        chunkEnd = min(chunkStart+CHUNKSIZE-1,nSampRead);
        
        fprintf('Reading samples %d:%d of %d...\n',firstsample+chunkStart-1,firstsample+chunkEnd-1,nSampRead);
        
        try
            data(:,chunkStart:chunkEnd) = read_mff_data(filename,'sample',firstsample+chunkStart-1,...
                firstsample+chunkEnd-1,1:mffhdr.nChans,dtype);
        catch err
            if strcmp(err.identifier,'MATLAB:Java:GenericException') && ...
                    ~isempty(strfind(err.message,'java.lang.OutOfMemoryError'))
                error(['%s\n\nreadegimff: Increase the size of your MATLAB Java heap space!\n' ...
                    'Find out how at http://www.mathworks.co.uk/support/solutions/en/data/1-18I2C/'],...
                    err.message);
            else
                rethrow(err);
            end
        end
        
        chunkStart = chunkEnd+1;
    end
    
    if ~isempty(evt)
        evt = evt(cell2mat({evt.sample}) >= firstsample & cell2mat({evt.sample}) <= lastsample);
        for e = 1:length(evt)
            evt(e).sample = evt(e).sample - firstsample + 1;
        end
    end
else
    try
        data = read_mff_data(filename,'epoch',1,mffhdr.nTrials,1:mffhdr.nChans,dtype);
    catch err
        if strcmp(err.identifier,'MATLAB:Java:GenericException') && ...
                ~isempty(strfind(err.message,'java.lang.OutOfMemoryError'))
            error(['%s\n\nreadegimff: Increase the size of your MATLAB Java heap space!\n' ...
                'Find out how at http://www.mathworks.co.uk/support/solutions/en/data/1-18I2C/'],...
                err.message);
        else
            rethrow(err);
        end
    end
    nSampRead = mffhdr.nSamples;
end

head.samp_rate = mffhdr.Fs;
head.nchan = mffhdr.nChans;
head.samples = nSampRead;
head.segments = mffhdr.nTrials;
head.segsamps = nSampRead;

javarmpath(mffjarpath);
