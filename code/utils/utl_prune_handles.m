function x = utl_prune_handles(x)
% Prune unreferenced workspace variables from anonymous functions in the given argument.
%
% In:
%   Data : an arbitrary data structure
%
% Out:
%   Data : the data structure with hidden references in anonymous functions removed
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-23

if iscell(x)
    % recurse into cell arrays...
    for k=1:numel(x)
        x{k} = utl_prune_handles(x{k}); end
elseif isstruct(x)
    % recurse into structs...
    for f=fieldnames(x)'
        if ~isempty(x)
            [x.(f{1})] = celldeal(utl_prune_handles({x.(f{1})})); end
    end
elseif isa(x,'function_handle')
    % inspect function handles
    f = char(x);
    if f(1) == '@'
        % anonymous function: take apart...
        parts = functions(x);
        % ... and put back together
        x = make_function(f,parts.workspace{1});
    end
end

% create a function handle
function f = make_function(decl__,workspace__)
% create workspace
for fn__=fieldnames(workspace__)'
    eval([fn__{1} ' = workspace__.(fn__{1}) ;']); end
clear workspace__ fn__;
% evaluate declaration
f = eval(decl__);
