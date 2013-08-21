% Sets the icaweights, icawinv fields under EEG structure to the
% chosen AMICA models'.
%
% Inputs:
% EEG      - EEG data structure
% model    - Model index of AMICA model to be chosen.
function [EEG, command] = pop_changeweights(EEG, model); 

    command = '';
    if nargin < 1
        help pop_loadmodout;
        return;
    end;
    
    command = '';
   if ~isfield(EEG.etc,'amica') || isempty(EEG.etc.amica)
    error('No AMICA solution found. You should first load AMICA components');
   end
   
   if nargin < 2
       
       allmodels= 1:EEG.etc.amica.num_models;
       
       if ~isfield(EEG.etc.amica,'modnames')
           EEG.etc.amica.modnames{length(allmodels)} = '';
       end
       
       str = [num2str(allmodels(1)) fastif(~isempty(EEG.etc.amica.modnames{1}) ,[' - ' EEG.etc.amica.modnames{1}],'')];
       for i = 2:EEG.etc.amica.num_models
           str = [str '|' num2str(allmodels(i)) fastif(~isempty(EEG.etc.amica.modnames{i}) ,[' - ' EEG.etc.amica.modnames{i}],'')];
       end
       editwhat2select = str;
       
       cb_enable = ['set(findobj(''parent'', gcbf, ''tag'', ''labelstring'')    , ''enable'', ''on'');'];
       cb_disable = ['set(findobj(''parent'', gcbf, ''tag'', ''labelstring'')    , ''enable'', ''off'');'];
       cb_label = ['if get(gcbo,''value''),' cb_enable ...
        'else,' cb_disable ...
        'end;'];
       
       uilist = {{'style' 'text' 'string' 'Choose model'}...
             {'style' 'popupmenu' 'string' editwhat2select 'value' 1}...
             {'style' 'checkbox' 'tag' 'labelcheck' 'string' 'Label the chosen model as:' 'value' 0 'callback' cb_label} ...
             {'style' 'edit' 'string' '' 'enable' 'off' 'tag' 'labelstring'}};
       uigeom = {[1 1] [1 1]};
       guititle = 'AMICA: pop_changeweights()';
       result = inputgui(uigeom, uilist, 'pophelp(''pop_modprobplotl'')', guititle, [], 'normal');
       if isempty(result)
            return;
       end
       model =result{1};
    end;
    
    if result{2}
          
           EEG.etc.amica.modnames{model} = result{3};
    end
    EEG.icaweights = EEG.etc.amica.W(:,:,model);  
    EEG.icasphere  = EEG.etc.amica.S(:,:); 
    EEG.icawinv    = EEG.etc.amica.A(:,:,model);
    EEG.icaact     = [];
    command = sprintf('EEG = pop_changeweights(%s,%d)',inputname(1),model);
    
    
    
