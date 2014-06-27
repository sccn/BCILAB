function res = utl_aggregate_results(varargin)
% Internal. Aggregate the given results (in any format allowed for ml_predict) into a single array.
% In:
%   Results... : results as produced by ml_predict
%
% Out:
%   Aggregate : all results concatenated into a single array
%
% See also:
%   ml_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-07

if isnumeric(varargin{1})
    try
        res = vertcat(varargin{:});
    catch e
        sizes = cellfun(@size,varargin,'UniformOutput',false);
        error('Failed to concatenate numeric prediction arrays -- make sure that the dimensionality of your predictions (%s) matches in each array; error message: %s.',hlp_tostring(sizes),e.message);
    end
elseif iscell(varargin{1}) && length(varargin{1}) > 1
    % make sure that the result formats are all equal
    v1 = varargin{1}([1 3:end]);
    for i=2:length(varargin)
        if ~isequal(v1,varargin{i}([1 3:end]))
            error('The formats of all prediction arrays to concatenate need to be identical.'); end
    end
    % grab the actual parameters
    params = cellfun(@(x)x{2},varargin,'UniformOutput',false);
    if isnumeric(varargin{1}{2})
        % numeric parameters are vertically concatenated
        try
            res = {varargin{1}{1} vertcat(params{:}) varargin{1}{3:end}};
        catch e
            sizes = cellfun(@size,params,'UniformOutput',false);
            error('Failed to concatenate prediction arrays -- make sure that the dimensionality of your distributions (%s) matches in each array; error message: %s.',hlp_tostring(sizes),e.message);
        end
    elseif iscell(varargin{1}{2})
        % cell parameters are horizontally concatenated
        try
            res = {varargin{1}{1} horzcat(params{:}) varargin{1}{3:end}};
        catch e
            sizes = cellfun(@size,params,'UniformOutput',false);
            error('Failed to concatenate prediction arrays -- make sure that the dimensionality of your distributions (%s) matches in each array; error message: %s.',hlp_tostring(sizes),e.message);
        end
    else
        error('This function can only concatenate prediction arrays of the form {X,Y, ...} where Y is either numeric or a cell array.');
    end
else
    error('This function can only concatenate arrays of predictions that are either all numeric or in the cell array form described in ml_predict.');
end

