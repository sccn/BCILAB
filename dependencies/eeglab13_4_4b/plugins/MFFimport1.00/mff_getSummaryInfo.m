function summaryInfo = mff_getSummaryInfo(filePath,datType)

if nargin <= 1
    datType = 0;
end

% create the MFFFile object
mfffileObj = javaObject('com.egi.services.mff.api.MFFFile', filePath, true);
if mfffileObj.isMFFFileLoaded();
    if datType == 0 %EEG
        [binObj blocks] = getEEGBlocks(mfffileObj, filePath);
    elseif datType == 1 %PIB
        [binObj blocks] = getPIBBlocks(mfffileObj, filePath);
    end
    numblocks = binObj.getNumberOfBlocks();
    if numblocks > 0 % Should always be
        blockObj = blocks.get(0); %zero based
        % Assumes sampling rate is the same in all blocks. Unpredictable if
        % does change. 
        sampRate = double(blockObj.signalFrequency(1)); % 1 based
        nChans = blockObj.numberOfSignals;
        [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0] = getEpochInfos(mfffileObj, filePath, sampRate);
        summaryInfo.sampRate = sampRate;
        summaryInfo.nChans = nChans;
        summaryInfo.binObj = binObj;
        summaryInfo.blocks = blocks;
        summaryInfo.epochType = epochType;
        summaryInfo.epochBeginSamps = epochBeginSamps;
        summaryInfo.epochNumSamps = epochNumSamps;
        summaryInfo.epochFirstBlocks = epochFirstBlocks;
        summaryInfo.epochLastBlocks = epochLastBlocks;
        summaryInfo.epochLabels = epochLabels;
        summaryInfo.epochTime0 = epochTime0;
        summaryInfo.mfffileObj = mfffileObj;
    else
        % Error: Signal has 0 blocks
    end
else
    % Error: Couldn't load MFF object
end
numblocks = binObj.getNumberOfBlocks();
for x = 0:numblocks-1
    blockObj = blocks.get(x);
    sampRate = double(blockObj.signalFrequency(1)); % 1 based
    nChans = blockObj.numberOfSignals;
    if sampRate ~= summaryInfo.sampRate;
        % Error: Inconsistent sampling rate
    end
    if nChans ~= summaryInfo.nChans;
        % Error: Inconsistent number of channels
    end
end

function [binObj blocks] = getEEGBlocks(mfffileObj, filePath)
% create Signal object and read in signal1.bin file
% !!!Replace with new call that returns EEG bin file. 
binfiles = mfffileObj.getBinFiles();
EEGBinInd = 0; 
EEGFiles = binfiles.elementAt(EEGBinInd);
binObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, EEGFiles, filePath);

% inspect Signal file.  How many blocks, samples, etc.
% this can be done without reading in any of the actual data
blocks = binObj.getSignalBlocks();

function [binObj blocks] = getPIBBlocks(mfffileObj, filePath)
% create Signal object and read in signal1.bin file
% !!!Replace with new call that returns PIB bin file. 
binfiles = mfffileObj.getBinFiles();
PIBBinInd = 1; 
PIBFiles = binfiles.elementAt(PIBBinInd);
binObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, PIBFiles, filePath);

% inspect Signal file.  How many blocks, samples, etc.
% this can be done without reading in any of the actual data
blocks = binObj.getSignalBlocks();

function [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0] = getEpochInfos(mfffileObj, filePath, sampRate)
epochList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs, 'epochs.xml', filePath);
epochListArray = epochList.getEpochs;
numEpochs = epochListArray.size;
epochBeginSamps = zeros(1,numEpochs);
epochNumSamps = zeros(1,numEpochs);
epochFirstBlocks = zeros(1,numEpochs);
epochLastBlocks = zeros(1,numEpochs);
epochNumSamps = zeros(1,numEpochs);
epochTime0 = zeros(1,numEpochs);
epochLabels = cell(1,numEpochs);
for p = 0:numEpochs-1
    anEpoch = epochListArray.get(p);
    epochBegin = anEpoch.getBeginTime();
    epochEnd = anEpoch.getEndTime();
    % Note: end time is first sample NOT in epoch. 
    epochBeginSamps(p+1) = mff_nanos2Sample(epochBegin, sampRate);
    epochTime0(p+1) = epochBeginSamps(p+1);
    epochNumSamps(p+1) = mff_nanos2Sample(epochEnd, sampRate) - epochBeginSamps(p+1);
    epochFirstBlocks(p+1) = anEpoch.getFirstBlock;
    epochLastBlocks(p+1) = anEpoch.getLastBlock;
    epochLabels{p+1} = 'epoch';
%fprintf('epoch diff %f\n', epochEnd - epochBegin)
end

% assumes the following, which should be true: 1-1 mapping between segments
% and epochs, including quantity and begin times.
% !! add checks and error cases? 
epochType = 'cnt';
categsFilesList = mfffileObj.getCategoriesFiles();
numCategsFiles = categsFilesList.size;
totalNumSegs = 0;
if numCategsFiles > 0 % should be either zero or one?
    categList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories, 'categories.xml', filePath);

    categListArray = categList.getCategories;
    numCategs = categListArray.size;
    for p = 0:numCategs-1
        aCateg = categListArray.get(p);
        categLabel = aCateg.getName;
        segList = aCateg.getSegments;
        numSegs = segList.size;
        totalNumSegs = totalNumSegs + numSegs;
        for q = 0:numSegs-1
            aSeg = segList.get(q);
            segBegin = aSeg.getBeginTime;
%segEnd = aSeg.getEndTime;
%fprintf('seg diff %f\n', segEnd - segBegin)
            segBeginSamp = mff_nanos2Sample(segBegin, sampRate);
            segInd = find(epochBeginSamps == segBeginSamp);
            epochLabels{segInd} = char(categLabel);
            time0 = aSeg.getEventBegin;
            time0Samp = mff_nanos2Sample(time0, sampRate);
            time0Samp = (time0Samp - segBeginSamp) + 1;
            epochTime0(segInd) = time0Samp;
        end
    end
    epochType = 'seg';
    % if epoch lengths are different, yet there are categories, than it's
    % var(iable) epoch type. 
    if size(unique(epochNumSamps),2) ~= 1 || size(unique(epochTime0),2) ~= 1
        epochType = 'var';
    end
end
%fprintf('totalNumSegs, numEpochs %d %d\n', totalNumSegs, numEpochs);
