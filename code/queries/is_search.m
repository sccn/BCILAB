function x = is_search(x)
% Determine whether some object is a search() expression.
x = all(isfield(x,{'head','parts'})) && numel(x)==1 && strcmp(char(x.head),'search');
