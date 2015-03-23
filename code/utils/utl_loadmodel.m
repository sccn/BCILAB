function model = utl_loadmodel(model)
% Load a given model.
% Model = utl_loadmodel(Model)
%
% In:
%   Model : The model to load; this can be a model struct, or a base workspace variable name, or a
%           file name, or a cell array of {file name, variable name} to refer to a variable inside a
%           file. Models are calibrated via bci_train or the GUI.
%
% Out:
%   Model : the loaded model structure
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2013-11-23

% parse the model argument
if ischar(model)
    % try to evaluate in the workspace
    try
        model = evalin('base',model);
    catch %#ok<CTCH>
        % if that fails, try to load it as a file name
        try
            model = io_load(model);
            if ~isfield(model,'tracking') || ~isfield(model.tracking,'prediction_function')
                % the loaded model is lacking the appropriate fields; check if there are variables
                % in the loaded data which are valid models
                candidates = {};
                for f = fieldnames(model)'
                    fname = f{1};
                    if isfield(model.(fname),'tracking') && isfield(model.(fname).tracking,'prediction_function')
                        candidates{end+1} = fname; end %#ok<AGROW>
                end
                if length(candidates) > 1
                    error('BCILAB:utl_loadmodel:ambiguous',['The file given as the model contains multiple candiate variables:\n' ...
                        hlp_tostring(candidates) '; please pass a file or model structure which is non-ambiguous.']); 
                elseif isempty(candidates)
                    error('BCILAB:utl_loadmodel:load_error','The given file contains no valid model.');
                else
                    model = model.(candidates{1});
                end
            end
        catch %#ok<CTCH>
            error('BCILAB:utl_loadmodel:load_error','The given model string could not be interpreted (neither as a file name nor as a workspace variable).');
        end
    end
elseif iscell(model) && length(model) == 2 && iscellstr(model)
    % two-element cell-string arrays are interpreted as {filename varname}.
    try
        model = getfield(io_load(model{1},model{2}),model{2});
    catch %#ok<CTCH>
        error('BCILAB:utl_loadmodel:load_error',['The file ' model{1} ' and/or its variable ' model{2} ' could not be loaded.']);
    end
elseif ~isstruct(model) || isempty(model)
    error('BCILAB:utl_loadmodel:invalid_model','The given data is not a valid model.');
elseif ~isfield(model,'tracking') || ~isfield(model.tracking,'prediction_function') || ~isfield(model.tracking,'filter_graph')
    error('BCILAB:utl_loadmodel:invalid_model','The given data structure is not a valid model (lacking required fields).');
end
