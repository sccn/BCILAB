function res = hlp_varargin2struct(args, varargin)
% struct = hlp_varargin2struct(varargin, defaults)
% converts a list of name-value pairs into a struct with values assigned to names
%
% In:
%   varargin: cell array of alternating names and values
%             instead of any name-value pair, a struct can be given, which is substituted in-place into the argument list
%   defaults: optional name-value list of defaults; if a default value is '<mandatory>', an error
%             is raised when no value is specified in varargin
%             instead of any name-value pair, a struct can be specified, which is substituted in-place into the argument list
%
% Notes:
%   some parameters may have multiple alternative names, which shall be remapped to the 
%   standard name within opts; alternative names are given together with the defaults,
%   by specifying a cell array of names instead of the name, as in the following example:
%   ... ,{'standard_name','alt_name_x','alt_name_y'}, default_value, ...          
%
% Out: 
%   struct: structure with values assigned to names
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-05

res = struct();
% it is allowed to specify just a struct in place of the varargin clause... (for further convenience)
if isstruct(args)    
    args = {args}; end
% splice all substructs into the name-value lists
args = flatten_structs(args);
varargin = flatten_structs(varargin);

% go through the defaults (varargin) and apply them (also, build a remapping table for additional names along the way)
default_name_for = struct();
for i=1:2:length(varargin)
    if iscell(varargin{i})
        % alternative names were specified, remember them
        for j=2:length(varargin{i})
            default_name_for.(varargin{i}{j}) = varargin{i}{1}; end
        varargin{i} = varargin{i}{1};
    end    
    try
        res.(varargin{i}) = varargin{i+1};
    catch
        error(['invalid field name specified in defaults: ' exp_fullform(varargin{i})]);
    end
end

% go through the arguments and override
for i=1:2:length(args)    
    if isfield(default_name_for,args{i})
        % rewrite alternative names
        args{i} = default_name_for.(args{i}); end
    try
        res.(args{i}) = args{i+1};
    catch
        error(['invalid field name specified in args: ' exp_fullform(args{i})]);
    end
end

% check for missing but mandatory args
for fn = fieldnames(res)'
    if strcmp(res.(fn{1}), '<mandatory>')
        error(['the parameter ' fn{1} ' was unspecified but is mandatory']); end        
end          


% substitute any structs in place of a name-value pair into the name-value list
function outargs = flatten_structs(inargs)
outargs = {};
k=1;
while k <= length(inargs)
    if isstruct(inargs{k})
        for fn=fieldnames(inargs{k})'
            outargs([end+1 end+2]) = {fn{1} inargs{k}.(fn{1})}; end
        k = k+1;
    else
        outargs([end+1 end+2]) = inargs([k k+1]);
        k = k+2;
    end
end
if mod(length(outargs),2)
    error('a non-matching number of names and values was specified.'); end