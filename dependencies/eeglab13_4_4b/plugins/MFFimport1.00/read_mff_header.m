function header = read_mff_header(filePath,datType)
summaryInfo = mff_getSummaryInfo(filePath,datType);

header.Fs = summaryInfo.sampRate;
header.nChans = summaryInfo.nChans;
nSamplesPre = 0;
if strcmp(summaryInfo.epochType, 'seg')
    header.nSamples = summaryInfo.epochNumSamps(1);
    header.nTrials = size(summaryInfo.epochBeginSamps,2);
    % if Time0 is the same for all segments...
    if size(unique(summaryInfo.epochTime0),2) == 1
        nSamplesPre = summaryInfo.epochTime0(1);
    end
else
    header.nSamples = sum(summaryInfo.epochNumSamps);
    header.nTrials = 1;
    header.nEpochs = size(summaryInfo.epochBeginSamps,2);
    header.epochStart = summaryInfo.epochBeginSamps;
    header.epochSize = summaryInfo.epochNumSamps;
end
sensorLayoutObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout, 'sensorLayout.xml', filePath);
sensors = sensorLayoutObj.getSensors();
nChans = 0;
for p = 1:sensors.size
    sensorObj = sensors.get(p-1); % sensors 0 based
    sensorType = sensorObj.getType;
    if sensorType == 0 || sensorType == 1
        tmpLabel = sensorObj.getName;
        if strcmp(tmpLabel,'')
            tmpLabel = sprintf('E%d', sensorObj.getNumber);
        end
        labels{p} = tmpLabel;
        nChans = nChans + 1;
    end
end
if nChans ~= header.nChans
    %Error
end
header.labels = labels;
