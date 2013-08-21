function res = is_needing_search(argform,args)
% test whether some argument pack requires a search or not (according to the specified argument format)

if strcmp(argform,'direct')
    % a search is specified by multielement arguments
    res = prod(max(1,cellfun(@length,args))) > 1;
elseif strcmp(argform,'clauses')
    % a search is specified by (possibly nested) search clauses
    res = contains_search(args);
else
    error('unsupported argument form.');
end

% test whether the given data structure contains a search clause
function res = contains_search(x)
if has_canonical_representation(x) && isequal(x.head,@search)
    res = true; 
elseif iscell(x)
    res = any(cellfun(@contains_search,x));
elseif isstruct(x) && numel(x) == 1
    res = contains_search(struct2cell(x));
else
    res = false;
end