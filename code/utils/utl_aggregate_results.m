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
        res = vertcat(varargin{:});
    elseif iscell(varargin{1}) && length(varargin{1}) > 1
        % make sure that the result formats are all equal
        v1 = varargin{1}([1 3:end]);
        for i=2:length(varargin)
            if ~isequal(v1,varargin{i}([1 3:end]))
                error('result formats are not identical.'); end
        end
        % grab the actual parameters
        params = cellfun(@(x)x{2},varargin,'UniformOutput',false);
        if isnumeric(varargin{1}{2})
            % numeric parameters are vertically concatenated
            res = {varargin{1}{1} vertcat(params{:}) varargin{1}{3:end}};
        elseif iscell(varargin{1}{2})
            % cell parameters are horizontally concatenated
            res = {varargin{1}{1} horzcat(params{:}) varargin{1}{3:end}};
        else
            error('unsupported result format.');
        end
    else
        error('unrecognized result format.');
    end
end
