% remove duplicates in a cell array of 
% <name,value> pairs (nvps)
% if duplicates exist, last pair in list is kept
% -------------------------------------------
function nvps = hlp_remDupNVPs(nvps, verbose)

if nargin<2
    verbose = false;
end

[tmp indices] = unique_bc(nvps(1:2:end));
if verbose && length(tmp) ~= length(nvps)/2
    fprintf('Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
end;
nvps = nvps(sort(union_bc(indices*2-1, indices*2)));