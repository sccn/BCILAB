
function [EEG com] = pop_selectmodel(EEG, model, type);

if nargin < 1
    help pop_selectmodel;
    return;
end;
com = '';

if ~isfield(EEG.etc,'amica') || isempty(EEG.etc.amica)
    error('No AMICA solution found. You should first load AMICA components');
end

guititle = 'Select/change data using AMICA model probabilities -- pop_selectmodel()';
data_options = 'Discard data for other models|Weight data by model probability';
if nargin < 2
    
    %newEEG = EEG;
    
    if EEG.etc.amica.num_models == 1
        error('Only one model exists');
    end
    
    % create user interface that allows choosing the model, weighting
    % method, selecting/discarding data. 
    allmodels= 1:EEG.etc.amica.num_models;
    
    if ~isfield(EEG.etc.amica,'modnames')
        EEG.etc.amica.modnames{length(allmodels)} = '';
    end
    
    str = [num2str(allmodels(1)) fastif(~isempty(EEG.etc.amica.modnames{1}) ,[' - ' EEG.etc.amica.modnames{1}],'')];
    for i = 2:EEG.etc.amica.num_models
        str = [str '|' num2str(allmodels(i)) fastif(~isempty(EEG.etc.amica.modnames{i}) ,[' - ' EEG.etc.amica.modnames{i}],'')];
    end
    editwhat2select = str;
    
    
    
    %smooth_options = 'None (default)|Hann';
    cb_enable = ['set(findobj(''parent'', gcbf, ''tag'', ''labelstring'')    , ''enable'', ''on'');'];
    cb_disable = ['set(findobj(''parent'', gcbf, ''tag'', ''labelstring'')    , ''enable'', ''off'');'];
    cb_label = ['if get(gcbo,''value''),' cb_enable ...
        'else,' cb_disable ...
        'end;'];
    
    uilist = { {'style' 'text' 'string' 'Select model'} ...
        {'style' 'popupmenu' 'string' editwhat2select 'value' 1} ...
        {'style' 'text' 'string' 'Data selection/weighting'} ...
        {'style' 'popupmenu' 'string' data_options 'value' 1} ...
        {'style' 'checkbox' 'tag' 'labelcheck' 'string' 'Label the chosen model as:' 'value' 0 'callback' cb_label} ...
        {'style' 'edit' 'string' '' 'enable' 'off' 'tag' 'labelstring'}};
    
    uigeom = {[1 2] [1 2] [1 1]};
    
    result = inputgui(uigeom, uilist, 'pophelp(''pop_selectmodel'')', guititle, [], 'normal');
    
    if isempty(result)
        
        return
    end
    
    
    
    
    model = result{1};
    
    
    if result{3}
        if ~isfield(EEG.etc.amica,'modnames')
            for i = 1:EEG.etc.amica.num_models
                EEG.etc.amica.modnames{i} = '';
            end
            
        end
        EEG.etc.amica.modnames{model} = result{4};
    end
    
    
    type = result{2};
    EEG= selectmodel(EEG,model,type);
    
else
    if ischar(model) && isfield(EEG.etc.amica,'mod_names')
        for i = 1:EEG.etc.amica.num_models
            if strcmpi(model,EEG.etc.amica.mod_names{i})
                model = i;
                break;
            end
            model = 0;
        end
    end
    if model < 1 
        error('Model not found');
    end    
    if model>=EEG.etc.amica.num_models
        error('Model number exceeds the total number of models');
    end
       
    
    if nargin<3
        uilist = {{'style' 'text' 'string' 'Data selection/weighting'}...
                 {'style' 'popupmenu' 'string' data_options 'value' 1}};
        uigeom = {[1 1.5]};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_selectmodel'')', guititle, [], 'normal');
        if isempty(result)
            return;
        end
        type = result{1};
    end
    EEG = selectmodel(EEG,model,type);
end


com = sprintf('%s = pop_selectmodel(%s,%d,%d)', inputname(1), inputname(1), model, type);
