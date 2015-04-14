function data = read_mff_data(filePath, indType, beginInd, endInd, chanInds, datType)
summaryInfo = mff_getSummaryInfo(filePath,datType);

if strcmp(indType, 'sample')
    [beginEpoch beginSample] = epochSample2EpochAndSample(beginInd, summaryInfo.epochNumSamps);
    [endEpoch endSample] = epochSample2EpochAndSample(endInd, summaryInfo.epochNumSamps);
    beginBlock = summaryInfo.epochFirstBlocks(beginEpoch);
    endBlock = summaryInfo.epochLastBlocks(endEpoch);
else
    beginBlock = summaryInfo.epochFirstBlocks(beginInd);
    endBlock = summaryInfo.epochLastBlocks(endInd);
end
data = read_mff_data_blocks(summaryInfo.binObj, summaryInfo.blocks, beginBlock, endBlock);
% if channel indeces were provided, downsample to the requested channels
if size(chanInds,1) ~= 0
    data = data(chanInds,:);
end
if strcmp(indType, 'sample')
    data = data(:,beginSample:beginSample + (endInd-beginInd));
elseif strcmp(summaryInfo.epochType, 'seg')
    nChans = size(data, 1);
    nSamples = summaryInfo.epochNumSamps(1);
    nTrials = (endInd - beginInd) + 1;
    data = reshape(data,nChans, nSamples, nTrials);
end

function data = read_mff_data_blocks(binObj, blocks, beginBlock, endBlock)
for blockInd = beginBlock-1:endBlock-1
    tmpdata = read_mff_data_block(binObj, blocks, blockInd);
    if blockInd == beginBlock-1
        data = tmpdata;
    else
        if size(data,1) == size(tmpdata,1)
            data = [data tmpdata];
        else
            % Error: blocks disagree on number of channels
        end
    end
end

function data = read_mff_data_block(binObj, blocks, blockInd)
blockObj = blocks.get(blockInd);
% to access the data for a block, it must be loaded first
blockObj = binObj.loadSignalBlockData(blockObj);
numChannels = blockObj.numberOfSignals;

% number of 4 byte floats is 1/4 the data block size
% That is divided by channel count to get data for each channel:
samplesTimesChannels = blockObj.dataBlockSize/4;
numSamples = samplesTimesChannels / numChannels;

% get first block, returned as bytes. 
data = blockObj.data;
% convert bytes to equivalent floating point values
data = typecast(data,'single');
data = reshape(data, numSamples, numChannels)';

function [epochNum sample] = epochSample2EpochAndSample(sampleNum, epochNumSamps)
epochNum = 1;
numSamps = epochNumSamps(epochNum);
while sampleNum > numSamps
    epochNum = epochNum + 1;
    numSamps = numSamps + epochNumSamps(epochNum);
end
numSamps = numSamps - epochNumSamps(epochNum);
sample = sampleNum - numSamps;
