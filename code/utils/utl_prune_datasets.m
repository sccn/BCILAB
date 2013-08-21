function x = utl_prune_datasets(x)
% Prune datasets from a given data structure.
%
% In:
%   Data : an arbitrary data structure
%
% Out:
%   Data : the data structure with data sets replaced by their expressions
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-12-06

if iscell(x)
    % recurse into cell arrays...
    for k=1:numel(x)
        x{k} = utl_prune_datasets(x{k}); end
elseif isstruct(x)    
    if all(isfield(x,{'data','chanlocs','srate'}))
        % this is a data set (or array thereof)
        if ~isscalar(x)
            % if array, process each one separately and combine results
            tmp = {};
            for k=1:numel(x)
                tmp{k} = utl_prune_datasets(x(k)); end
            x = reshape([tmp{:}],size(x));
        else
            % if scalar
            if isfield(x,'tracking') && isfield(x.tracking,'expression')
                % replace by its expression
                x = x.tracking.expression;
            else
                % or get rid of it entirely
                x = '<EEGLAB dataset: pruned>';
            end
        end
    else
        % recurse into structs...
        for f=fieldnames(x)'
            if ~isempty(x)
                [x.(f{1})] = celldeal(utl_prune_datasets({x.(f{1})})); end
        end
    end
end
