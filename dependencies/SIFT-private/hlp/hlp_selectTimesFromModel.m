function MODEL = hlp_selectTimesFromModel(MODEL,timesInSeconds,mode)
% select a range of time indices from a VAR Model (relative to
% MODEL.winStartTimes
%
% Tim Mullen, SCCN/INC UCSD 2012

if nargin<3
    mode = 'sequence';
end

timeIndices = getindex(MODEL.winStartTimes,timesInSeconds);

if length(timeIndices)==2 || strcmpi(mode,'range')
    timeIndices = timeIndices(1):timeIndices(end);
end

fnames = fieldnames(MODEL);
for f = 1:length(fnames)
    if ndims(MODEL.(fnames{f}))==4
        % prune the third dimension (time)
%         MODEL.(fnames{f}) = MODEL.(fnames{f})(:,:,timeIndices,:);
    elseif iscell(MODEL.(fnames{f})) ...
        && length(MODEL.(fnames{f}))==length(MODEL.winStartTimes)
    
        % prune cell arrays
        MODEL.(fnames{f}) = MODEL.(fnames{f})(timeIndices);
    end
end

MODEL.winStartTimes = MODEL.winStartTimes(timeIndices);