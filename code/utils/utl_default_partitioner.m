function res = utl_default_partitioner(data,inds,varargin)
% The default partitioner for generic data (used in cross-validations). 
%   if only data and no inds are given, return the number of available indices
%   otherwise return data restricted to the given index set
%
% See also:
%   utl_crossval, utl_partition_bundle, set_partition

if ~exist('inds','var')
    inds = []; end

% check supported special cases first 
if iscell(data) && length(data) == 2 && isnumeric(data{1}) 
    % these are the std formats used by ml_train and ml_predict
    if isnumeric(data{2}) && size(data{1},1) == size(data{2},1) && size(data{1},2) ~= 1 && size(data{2},2) == 1
        % a cell array of the form {NxF,Nx1} -- indicating data and target variable; partitioned along N; F>1
        if isempty(inds)
            res = length(data{2});
        else
            res = {data{1}(inds,:),data{2}(inds,:)};
        end
    elseif iscell(data{2}) && length(data{2}) > 1 && ischar(data{2}{1}) && ((isnumeric(data{2}{2}) && size(data{1},1) == size(data{2}{2},1)) || (iscell(data{2}{2}) && size(data{1},1) == length(data{2}{2})))
        % a cell array of the form {NxF,{distrib,[NxD],...}} or {NxF,{distrib,{N},...}} --  indicating data and target variable, partitioned along N
        if isempty(inds)
            res = size(data{1},1);
        else
            if isnumeric(data{2}{2})
                res = {data{1}(inds,:), {data{2}{1},data{2}{2}(inds,:), data{2}{3:end}}};
            else
                res = {data{1}(inds,:), {data{2}{1},data{2}{2}(inds), data{2}{3:end}}};
            end
        end
    end
elseif isstruct(data) && isscalar(data) && isfield(data,'X')   
    % a struct with fields X & ...
    if isfield(data,'Y') && size(data.X,1) == size(data.Y,1)
        % ... upper-case Y: partition X and Y
        if isempty(inds)
            res = ult_default_partitioner(data.Y);
        else
            res = data;
            res.X = ult_default_partitioner(res.X,inds);
            res.Y = ult_default_partitioner(res.Y,inds);
        end
    elseif isfield(data,'y') && size(data.X,1) == size(data.y,1)
        % ... lower-case y: partition X and y
        if isempty(inds)
            res = ult_default_partitioner(data.y);
        else
            res = data;
            res.X = ult_default_partitioner(res.X,inds);
            res.y = ult_default_partitioner(res.y,inds);
        end
    end
elseif isstruct(data) && isfield(data,{'data','chanlocs'})
    % EEGLAB data set
    if isempty(inds)
        res = exp_eval(set_partition(data,[],varargin{:}),inf);
    else
        res = set_partition(data,inds,varargin{:});
    end
elseif isstruct(data) && isfield(data,{'streams'})
    % stream bundle
    res = utl_partition_bundle(data,inds,varargin{:});
end


if ~exist('res','var')
    % no special case: get the highest non-singleton dimension in the data
    dim = find(size(data)~=1,1,'last');
    if isnumeric(data)
        % in the numeric array case, we never partition along the first dimension
        dim = max(2,dim); end
    
    if isempty(inds) 
        res = size(data,dim);
    else
        switch dim
            case 1
                res = data(inds);
            case 2
                res = data(:,inds);                
            case 3
                res = data(:,:,inds);
            case 4
                res = data(:,:,:,inds);
            otherwise
                error('More than 4-dimensional data is not (yet) supported.');
        end
    end
end
