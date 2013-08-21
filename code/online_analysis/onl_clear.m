function onl_clear()
% Clear all online streams and predictors.
%
% Multiple online streams can be running in the background and consume resources. This function 
% allows to clear the respective resources (streams and predictors).
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-18

% find all variables
vars = evalin('base','whos');
vars = {vars(strcmp({vars.class},'struct')).name};
to_delete = {};
% for all struct vars...
for v=vars
    % select those that are streams or predictors
    var = evalin('base',v{1});
    if isfield(var,{'buffer','smax'})
        to_delete{end+1} = v{1}; end
    if all(isfield(var,{'pipelines','tracking'})) && isfield(var.tracking,'prediction_function')
        to_delete{end+1} = v{1}; end    
end

% delete them
if ~isempty(to_delete)
    evalin('base',['clear ' sprintf('%s ',to_delete{:}) ';']); end
