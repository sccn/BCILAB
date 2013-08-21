function EEG = pop_runamica(EEG)
amicaVersions = {'AMICA10'};
computeVersions = {'In parallel' 'Locally'};
shareCompOptions = {'No' 'Yes'};

callbackAMICAVersion = ['if get(gcbo,''value'') ~= 1,' ...
    'set(findobj(''parent'',gcbf,''tag'',''shareComp''),''enable'',''off'');'...
    'else,' ...
    'set(findobj(''parent'',gcbf,''tag'',''shareComp''),''enable'',''on'');'...' ...
    'end;'];

defaultOutputDirectory = [pwd '/amicaout'];
uilist = {{'style' 'text' 'string' 'AMICA version'} ...
    {'style' 'popupmenu' 'string' amicaVersions 'value' 1 'callback' callbackAMICAVersion} ...
    {'style' 'text' 'string' 'Output directory'} ...
    {'style' 'edit' 'string' defaultOutputDirectory} ...
    {'style' 'text' 'string' 'Run AMICA'} ...
    {'style' 'popupmenu' 'string' computeVersions 'value' 1} ...
    {'style' 'text' 'string' 'Number of Models'} ...
    {'style' 'edit' 'string' '1'} ...
    {'style' 'text' 'string' 'Number of Processors'} ...
    {'style' 'edit' 'string' '4'} ...
    {'style' 'text' 'string' 'Max. number of iterations'} ...
    {'style' 'edit' 'string' '2000'} ...
    {'style' 'text' 'string' 'Share components across models'} ...
    {'style' 'popupmenu' 'string' shareCompOptions 'value' 1 'tag' 'shareComp' 'enable' 'on'} ...
    {'style' 'text' 'string' 'Additional commandline options'} ...
    {'style' 'edit' 'string' ''}};

uigeom = {[1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1]};
guititle = 'Run AMICA -- pop_runamica()';
result = inputgui(uigeom, uilist, 'pophelp(''pop_runamica'')', guititle, [], 'normal');
GUIOptionKeywords = {'outdir' 'qsub' 'num_models' 'numprocs' 'max_iter' 'share_comps'};
if isempty(result)
    return
end

for i = 1:length(GUIOptionKeywords)
    keyword = GUIOptionKeywords{i};
    arglist{2*i-1} = keyword;
    if strcmpi(keyword,'outdir')
        arglist{2*i} = result{2}; 
    end
    if strcmpi(keyword,'qsub')
        arglist{2*i} = fastif(result{3}==1,'on','off');
    end
    if strcmpi(keyword,'num_models')
        arglist{2*i} = str2num(result{4});
    end
    if strcmpi(keyword,'numprocs')
        arglist{2*i} = str2num(result{5});
    end
    if strcmpi(keyword,'max_iter')
        arglist{2*i} = str2num(result{6});
    end
    if strcmpi(keyword,'share_comps')
        arglist{2*i} = fastif(result{7}==1,0,1);
    end
end

additionalOptions = eval( [ '{' result{8} '}' ]);
for i = 1:length(additionalOptions)
    arglist{end+1} = additionalOptions{i};
end

amicaVersion = result{1};

if amicaVersion == 1
    if EEG.trials>1
        tmpdata = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
        mod = runamica10(tmpdata,[],size(tmpdata,1),size(tmpdata,2),arglist{:});
    else
        mod = runamica10(EEG.data,[],size(EEG.data,1),size(EEG.data,2),arglist{:});
    end
    fprintf('AMICA output is going to be in the folder %s',outdir);
else
end



