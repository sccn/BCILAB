function result = combine_structs(args,structmask)
% Turn a list of name-value pairs/structs into a single struct.
% Struct = combine_structs(Arguments,StructMask)
%
% In:
%   Arguments : cell array of function arguments in form of name-value pairs, possibly interleaved
%               with structs {'name',value,STRUCT,'name',value,'name',value, ...}
%
%   StructMask : a bitmask that encodes at what positions Arguments contains structs; optional.
%
% Out:
%   Struct : A struct with values assigned to fields
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-03

if nargin < 2
    structmask = cellfun('isclass',args,'struct'); end

if isscalar(structmask)
    % inputs has already a single struct
    result = args{1};
else
    % obtain flat NVP list
    nvps = flatten_structs(args,structmask);
    % split names/values and find indices of the last assignment for each name
    nvps = reshape(nvps,2,[]);
    [s,indices] = sort(nvps(1,:));
    indices(strcmp(s(1:end-1),s(2:end))) = [];
    % build the struct
    result = cell2struct(nvps(2,indices),nvps(1,indices),2);
end
