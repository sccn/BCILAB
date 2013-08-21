function y = utl_get_argument(x,name)
% Get an argument (by name) of an expression that involves name-value pairs.
% Value = utl_get_argument(Expression,Name)
%
% In:
%   Expression : an expression with name-value arguments, for example a node of a 
%                filter graph
%
%   Name : name of the argument to extract
%
% Out:
%   Value : value of the desired argument
%
% See also:
%   utl_find_filter
%
% Examples:
%   % first get a filter node from a model
%   ica_node = utl_find_filter(lastmodel,'flt_ica')
%   % then extract its 'variant' argument, which is the ICA variant used
%   variant = utl_get_argument(ica_node,'variant')
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2012-01-12

if ~isfield(x,'parts')
    error('The given object is not an expression.'); end
match = find(strcmp(x.parts,name),1,'last');
if isempty(match)
    if ~iscellstr(x.parts(1:2:end))
        error('The given expression does not have name-value pairs as arguments.');
    else
        error('The given name is not an argument of the function.');
    end
end
y = x.parts{match+1};
