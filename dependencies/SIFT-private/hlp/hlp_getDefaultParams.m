% returns the default parameters in structure params for a given prefix
% class of functions
% see hlp_getDefaultArglist.m
function params = hlp_getDefaultParams(prefix)

params = hlp_getDefaultArglist(prefix);
params = params(:,[1 4])';
params = hlp_varargin2struct(params(:));
