% determine whether some object is a raw EEGLAB data set with no BCILAB constituents
function res = is_raw_dataset(x)
res = all(isfield(x,{'data','srate'})) && ~isfield(x,'tracking');
