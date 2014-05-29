function block = ns_read(h)
% Block reader for Neuroscan Scan Recorder
% ns_read(h)
%
%
% In:
%   h : handle to an existing Neuroscan connection
%
% Out:
%   block: cell array data block containing EEG and event values
%     
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald

data = [];
if (~h.initialized)
    return;
end

try
    % check for existing data in socket buffer
    headerBytes = pnet(h.handle, 'read', h.headerSize, 'uint8', 'noblock');
    
    events = [];
                
    if ~isempty(headerBytes)
        
        headerData = ns_parseheader(headerBytes);
        
        if strcmp(headerData.id, 'DATA')
            
            % wait until the data packet body is available
            dataPreview = [];
            while length(dataPreview) < headerData.bodysize
                % wait until the rest of the packet is available
                dataPreview = pnet(h.handle, 'read', headerData.bodysize,...
                    'uint8', 'view', 'noblock');
            end
            
            % read the data packet body
            dataBytes = pnet(h.handle, 'read', headerData.bodysize, 'uint8', 'noblock');
            
            dataBytes = reshape(uint8(dataBytes), [h.bytesPerSample, h.totalChan, h.samplesPerBlock]);
            
            markerValues = squeeze((dataBytes(1,h.markerChanIdx, :)));
            markerPoints = find(markerValues);
                        
            if ~isempty(markerPoints)
                markers = markerValues(markerPoints);
                latencies = markerPoints;
                for m = 1:length(markers)
                   events(m).type = num2str(markers(m)); %#ok<AGROW>
                   events(m).latency = latencies(m); %#ok<AGROW>
                end
            else
               events = []; 
            end
            
            dataBytes(:,h.markerChanIdx, :) = [];
            
            dataCell = squeeze(num2cell(dataBytes, 1));
            dataBlock = cellfun(@(x) typecast(x, h.datatype), dataCell);
            data = horzcat(data, dataBlock);
                        
        else
            disp('unknown message');
        end
                
    end
    
        
catch er
    disp(er.message);
    rethrow(er);
end

if ~isempty(data)
    % scale data to uV units
    data = bsxfun(@times,data,h.resolution);
end

block = {data,events};

