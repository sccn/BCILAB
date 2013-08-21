function com = pop_topohistplot2(EEG, components, model)
com = '';
if nargin < 1
    help pop_topohistplot
end
try empty = isempty(EEG.etc.amica);
    if empty
        error('No AMICA solution found. You should first load AMICA components');
    end
        if nargin < 2
            
            editwhat2plot = ['1:' int2str(size(EEG.icaweights,1))];
            allmodels= 1:EEG.etc.amica.num_models;
            editwhat2select(1) = num2str(allmodels(1));
            for i = 2:EEG.etc.amica.num_models
                editwhat2select(2*i-2) = '|';
                editwhat2select(2*i-1) = num2str(allmodels(i));
            end
            
            uilist = {{'style' 'text' 'string' 'Select components'}...
                {'style' 'edit' 'string' editwhat2plot}...
                {'style' 'text' 'string' 'Select model'} ...
                {'style' 'popupmenu' 'string' editwhat2select 'value' 1}...
                {}...
                {}}; 
            uigeom = {[1 1] [1 1]};
            guititle = 'Select components and model -- pop_topohistplot()';
            result = inputgui(uigeom, uilist, 'pophelp(''pop_topohistplot'')', guititle, [], 'normal');
            if isempty(result)
                return
            end
            components = eval( [ '[' result{1} ']' ] );
            model = result{2};
            if length(components) > EEG.nbchan
                tmpbut = questdlg2(...
                  ['This involves drawing ' int2str(length(components)) ' plots. Continue ?'], ...
                         '', 'Cancel', 'Yes', 'Yes');
            if strcmp(tmpbut, 'Cancel'), return; end;
            end;
            if isempty(components), error('Nothing to plot; enter parameter in first edit box'); end;
            
            
        else
            if nargin < 3
                model = 1;
            end
        end
        
        h = model;
        alpha = EEG.etc.amica.alpha;
        sbeta = EEG.etc.amica.sbeta;
        rho = EEG.etc.amica.rho;
        mu = EEG.etc.amica.mu;
        nbgraph = length(components);
        rowcols(2) = ceil(sqrt(nbgraph));
        rowcols(1) = ceil(nbgraph/rowcols(2));
        
        s = -10:0.01:10;
        for index = 1:length(components)
            sum = 0;
            for j = 1:size(alpha,1)
                sum = sum + alpha(j,components(index),h)*sbeta(j,components(index),h)*(1/(2*gamma(1+1/rho(j,components(index),h))))*exp(-abs(sbeta(j,components(index),h)*(s-mu(j,components(index),h)).^rho(j,components(index),h)));
            end
            p(index,:) = sum;
            subplot(rowcols(1),2*rowcols(2),2*index-1); topoplot(EEG.etc.amica.A(:,components(index),h),EEG.chanlocs);
            drawnow
            subplot(rowcols(1),2*rowcols(2),2*index); plot(s,p(index,:));
            drawnow;
        end
        
        
        
catch
    error('No AMICA solution found. You should first load AMICA components');
end
        
        
        
        
