function args = flatten_structs(args,structmask)
% Turn a list of name-value pairs/structs into a flat list of name-value pairs
% NVPs = flatten_structs(Arguments,StructMask)
%
% In:
%   Arguments : cell array of function arguments in form of name-value pairs, possibly interleaved
%               with structs {'name',value,STRUCT,'name',value,'name',value, ...}
%
%   StructMask : a bitmask that encodes at what positions Arguments contains structs; optional.
%
% Out:
%   NVPs : a cell array of name-value pairs, where the structs have been flattened into 
%          name-value pairs, too
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-11-03

if nargin < 2
    structmask = cellfun('isclass',args,'struct'); end

if any(structmask)
    persistent indexcache; %#ok<TLEV>
    try
        % try to look up splice points from cache
        if length(structmask) < 62
            field = char('a'+structmask);
        else
            str = java.lang.String(char('a'+structmask));
            field = ['h' num2str(str.hashCode()+2^31)];
        end
        splicepos = indexcache.(field);
    catch %#ok<CTCH>
        % pre-calculate splice points and cache results
        splicepos = [];
        k = 1;
        while k <= length(args)
            if isstruct(args{k})
                splicepos(end+1) = k; %#ok<AGROW>
                k = k+1; % struct case
            else
                k = k+2; % NVP case
            end
        end
        splicepos = splicepos(end:-1:1);
        indexcache.(field) = splicepos;
    end
    
    % splice structs in
    for k = splicepos
        args = [args(1:k-1) reshape([fieldnames(args{k}) struct2cell(args{k})]',1,[]) args(k+1:end)]; end
end
