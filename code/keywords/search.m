function res = search(varargin)
%   Special keyword to declare parameter values to search over.
%   search(x,y,z,...) or search([x,y,z,....]) or search({x,y,z,...}) expresses
%   a parameter range (i.e. set of possible parameter values, to be searched for) in a grid search.
%
% Note:
%   search() may also be nested inside a data structure, to indicate that all instances of that data structure, 
%   with each instances having one of the specified values substituted in it, shall be searched.
%
%   search() returns an expression data structure that encapsulates the specifies search range,
%   for use as argument to utl_gridsearch, utl_searchmodel, utl_nested_crossval, etc.
%
% See also:
%   utl_gridsearch, utl_searchmodel, utl_nested_crossval
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22


if length(varargin) == 1
    if iscell(varargin{1})
        varargin = varargin{1};
    elseif isnumeric(varargin{1})
        varargin = num2cell(varargin{1});
    end
end

res = struct('head',{@search},'parts',{varargin});