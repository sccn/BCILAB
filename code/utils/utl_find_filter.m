function [y,s] = utl_find_filter(x,name)
% Extract a filter node of a given name (e.g. flt_ica) from a model or filter graph.
% function Filter = utl_find_filter(Model,Name)
%
% This function is useful to get access to the filter parameters of a model, e.g. for visualiation.
%
% In:
%   Model : the model (or filter graph, or expression) in which to search
%   
%   Name : name of the filter stage to search (this is the function name)
%
% Out:
%   Filter : the extracted filter node (in particular its set of arguments)
%
%   Subscript : the subscript where the filter was found (usable with subsref or subsasgn)
%
% See also:
%   utl_get_argument
%
% Examples:
%   % get the ICA filter node of a model that used ICA
%   ica_node = utl_find_filter(lastmodel,'flt_ica')
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2012-01-12

%  S is a structure array with the fields:
%         type -- string containing '()', '{}', or '.' specifying the
%                 subscript type.
%         subs -- Cell array or string containing the actual subscripts.
 
s = [];
y = [];
if isfield(x,'tracking') && isfield(x.tracking,'filter_graph')
    % descend into models
    [y,s] = utl_find_filter(x.tracking.filter_graph,name);
    s = [struct('type',{'.','.'},'subs',{'tracking','filter_graph'}) s];
elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
    % descend into signals
    [y,s] = utl_find_filter(x.tracking.expression,name);
    s = [struct('type',{'.','.'},'subs',{'tracking','expression'}) s];
elseif iscell(x)
    % descend into cell array
    for k=1:length(x)
        [y,s] = utl_find_filter(x{k},name);
        if ~isempty(y)
            s = [struct('type',{'{}'},'subs',{{k}}) s]; %#ok<AGROW>
            break; 
        end
    end
elseif all(isfield(x,{'head','parts'}))
    if ~strcmp(char(x.head),name)
        % scan parts
        [y,s] = utl_find_filter(x.parts,name);
        s = [struct('type',{'.'},'subs',{'parts'}) s];
    else
        % found it
        y = x;
        s = struct('type',{'()'},'subs',{{1}});
    end
end
