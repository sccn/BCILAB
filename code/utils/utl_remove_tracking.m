function X = utl_remove_tracking(X)
% Remove the .tracking field (if present) from a data structure
%
% In:
%   X : the data structure
%
% Out:
%   X : data structure without the tracking field.
%
% Notes:
%   Removing this field clears any history of the data structure and so will be treated
%   by BCILAB as a "raw" EEGLAB data set. One case where this is useful is when a series of 
%   filters and dataset editing steps is applied manually to the data and then the data should be
%   processed as if it were raw (for example to shut off errors about non-causal processing steps
%   applied online). Obviously this is very rarely a good idea.

if isfield(X,'tracking')
    X = rmfield(X,'tracking'); end