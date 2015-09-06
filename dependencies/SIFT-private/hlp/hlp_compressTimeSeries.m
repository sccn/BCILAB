function signal = hlp_compressTimeSeries(signal,tsfields,newtype,remfields)
% tsfields: a cell array of fields to 'compress' by averaging across cols
%           (after taking the absolute value)
% signal: an eeg data structure
% newtype: new datatype to cast to
% remfields: fields to remove

if nargin<3
    newtype = '';
end
if nargin<4
    rmfields = {};
end

for fn=tsfields
    f = fn{1};
    signal.(f)     = abs(signal.(f));
    signal.(f)     = mean(signal.(f),2);
    if ~isempty(newtype)
        signal.(f) = cast(signal.(f),newtype);
    end
end

signal = rmfield(signal,remfields);

signal.xmax          = signal.xmin;
signal.pnts          = 1;
