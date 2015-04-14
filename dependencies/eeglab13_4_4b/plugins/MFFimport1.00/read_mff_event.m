function events = read_mff_event(filePath)

events = [];
summaryInfo = mff_getSummaryInfo(filePath);
epochBeginSamps = summaryInfo.epochBeginSamps;
epochNumSamps = summaryInfo.epochNumSamps;
epochFirstBlocks = summaryInfo.epochFirstBlocks;
epochLastBlocks = summaryInfo.epochLastBlocks;

infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);
beginTime = infoObj.getRecordTime();

% Get event file names. 
% Need to blend in events from different tracks
% Sort by sample num (only should be necessary w multiple tracks
eventtracknamelist = summaryInfo.mfffileObj.getEventTrackFiles();
eventtrackcount = eventtracknamelist.size();

if eventtrackcount > 0
    MFFUtil = javaObject('com.egi.services.mff.utility.MFFUtil');
    % read through each event track and inspect count of events
    % and beginTime for last event in each
% Code to count events before putting them into the array so that array can
% be pre-allocated. Didn't follow up, for now. !!Issues: 
% -- How to preallocate array of structures.
% -- Doubles the javaObject creates, which might have memory leak issues.
%     totalNumEvents = 0;
%     for tracknum = 0:eventtrackcount-1
%         eventTrackObj = mff_getObject('com.egi.services.mff.api.EventTrack', eventtracknamelist.get(tracknum), filePath);
%         eventList = eventTrackObj.getEvents;
%         totalNumEvents = totalNumEvents + eventList.size;
%     end
    eventInd = 0;
    for p = 1:size(summaryInfo.epochBeginSamps,2)
        eventInd = eventInd + 1;
        events(eventInd).type = ['break ' summaryInfo.epochType];
        events(eventInd).sample = samples2EpochSample(summaryInfo.epochBeginSamps(p), epochBeginSamps, epochNumSamps);
        events(eventInd).value = summaryInfo.epochLabels{p};
        events(eventInd).offset = [];
        events(eventInd).duration = summaryInfo.epochNumSamps(p);
        events(eventInd).timestamp = []; % or calculate string
        eventInds(eventInd,1) = events(eventInd).sample;
        eventInds(eventInd,2) = eventInd;
        if strcmp(summaryInfo.epochType, 'var')
            eventInd = eventInd + 1;
            events(eventInd).type = 't0';
            events(eventInd).sample = samples2EpochSample((summaryInfo.epochTime0(p) + summaryInfo.epochBeginSamps(p)) - 1, epochBeginSamps, epochNumSamps);
            events(eventInd).value = 't0';
            events(eventInd).offset = [];
            events(eventInd).duration = 1;
            events(eventInd).timestamp = []; % or calculate string
            eventInds(eventInd,1) = events(eventInd).sample;
            eventInds(eventInd,2) = eventInd;
        end
    end
    for tracknum = 0:eventtrackcount-1
        eventTrackObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack, eventtracknamelist.get(tracknum), filePath);
        eventList = eventTrackObj.getEvents;
        numEvents = eventList.size;
        for p = 0:numEvents-1
            % Java arrays are 0 based
            theEvent = eventList.get(p);
            eventTime = theEvent.getBeginTime;

            %Need to convert to samples (get sampling rate up above)
            eventTimeInNS = MFFUtil.getTimeDifferenceInNanoseconds(eventTime , beginTime);
            eventTimeInSamples = mff_nanos2Sample(eventTimeInNS, summaryInfo.sampRate);
            eventTimeInEpochSamples = samples2EpochSample(eventTimeInSamples, epochBeginSamps, epochNumSamps);
            if eventTimeInEpochSamples < 1
                % Error: invalid sample number
            else
                eventInd = eventInd + 1;
%                fprintf('%d %d %d\n', eventTimeInSamples, eventTimeInEpochSamples, eventTimeInSamples-eventTimeInEpochSamples);
                % Matlab arrays are 1 based
                events(eventInd).type = char(theEvent.getCode);
                events(eventInd).sample = eventTimeInEpochSamples;
                events(eventInd).value = [];
                events(eventInd).offset = [];
                events(eventInd).duration = nanos2SampleDuration(theEvent.getDuration, summaryInfo.sampRate);
                events(eventInd).timestamp = []; %eventTime;
                events(eventInd).codes = {};
                
                keyList = theEvent.getKeys;
                for key = 1:keyList.size
                    thisKey = keyList.elementAt(key-1);
                    events(eventInd).codes{key,1} = (char(thisKey.getCode));
                    events(eventInd).codes{key,2} = str2double(char(thisKey.getData));
                    if isnan(events(eventInd).codes{key,2})
                        events(eventInd).codes{key,2} = char(thisKey.getData);
                    end
                end
                
                eventInds(eventInd,1) = events(eventInd).sample;
                eventInds(eventInd,2) = eventInd;
            end
        end
    end
    eventInds = sortrows(eventInds);
    for p = 1:eventInd
        nextEventInd = eventInds(p,2);
        sortedEvents(p) = events(nextEventInd);
    end
    events = sortedEvents;
end

function sampleNum = nanos2SampleDuration(ns, sampRate)
sampDuration = 1000000000/sampRate;
sampleNum = ns/sampDuration;
%fprintf('%f\n', sampleNum);
sampleNum = ceil(sampleNum);
%fprintf('%d\n', sampleNum);

function epochSampleNum = samples2EpochSample(sampleNum, epochBeginSamps, epochNumSamps)
numEpochs = size(epochBeginSamps,2);
p = 1;
epochSampleNum = 0;
while (sampleNum > (epochBeginSamps(p) + epochNumSamps(p) - 1)) && (p < numEpochs)
    epochSampleNum = epochSampleNum + epochNumSamps(p);
    p = p+1;
end
if p <= numEpochs
    if sampleNum >= epochBeginSamps(p)
        epochSampleNum = epochSampleNum + ((sampleNum - epochBeginSamps(p))+1);
    else
        epochSampleNum = -1;
        % Error: sample falls between epochs
    end
else
    epochSampleNum = -2;
    % Error: sample is after last epoch
end



