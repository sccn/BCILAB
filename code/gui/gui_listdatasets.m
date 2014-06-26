function [list,flat,indexable] = gui_listdatasets()
% Compile a list of all datasets that are currently in the workspace
% List = gui_listdatasets()
%
% Out:
%   List : cell array of {'Groupname', {data1,data2,data3, ...}, 'Groupname', {data1,data2}, ...}
%
%   Flat : pretty-printed cell-string array of List contents
%  
%   Indexable : a flat cell array of datas ets that can be indexed in accordance with Flat
%
% Notes:
%   if no output is taken, the function displays the Flat list.
%
%   The following groups are generated:
%   'Loaded via BCILAB' --> expressions which contain an io_loadset node somewhere...
%   'Loaded via EEGLAB' --> full continuous (non-empty) datasets, including EEG
%
%   The datasets are variable names of the form 'varname(index) (comments)', where (index) 
%   is omitted for one-element struct arrays (i.e. regular structs). Thus, the command
%     evalin('base',strtok(name)) gives the appropriate workspace variable.
%
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

list = {};

% --- aggregate all EEGLAB datasets in the workspace ---
datasets = {};

vars = evalin('base','whos');
% go through all structs in the workspace
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    % check if it is an EEGLAB data set (array)
    if all(isfield(var,{'data','chanlocs','srate','xmin'})) 
        % go through all entries
        for k=1:length(var)
            % check if it is non-empty and continuous
            if ~isempty(var(k).data) && ndims(var(k).data) == 2
                % list it
                if length(var) == 1
                    varname = vars(v).name;
                else
                    varname = [vars(v).name '(' num2str(k) ')'];
                end
                if isfield(var(k),'setname') && ~isempty(var(k).setname)
                    name = [varname ' ("' var(k).setname '")'];
                else
                    name = varname;
                end
                datasets{end+1} = name;
            end
        end
    end
end

list = [list {'Loaded via EEGLAB', datasets}];


% --- aggregate all BCILAB datasets in the workspace ---
datasets = {};

vars = evalin('base','whos');
% go through all structs in the workspace
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if isfield(var,{'head','parts'})
        % found an expression: check if it contains an io_loadset() subexpression
        occurrences = utl_cases(var,exp_blank(@io_loadset),Inf);
        if ~isempty(occurrences)
            % list it
            if length(occurrences) == 1
                if isfield(occurrences{1},'tracking') && isfield(occurrences{1}.tracking,'expression')
                    occurrences{1} = occurrences{1}.tracking.expression; end
                % there is exactly 1 io_loadset(filepath, ...) in the expression: put the filename in brackets
                [path,file,ext] = fileparts(occurrences{1}.parts{1}); %#ok<ASGLU>
                if isequal(var.head,@io_loadset)
                    name = [vars(v).name ' ("' file ext '")']; 
                else
                    name = [vars(v).name ' ("' file ext '"; modified)']; 
                end
            else
                name = vars(v).name;
            end
            datasets{end+1} = name;
        end
    end
end

list = [list {'Loaded via BCILAB', datasets}];

[flat,indexable] = gui_print_grouplist(list);

% print it
if nargout == 0
    fprintf('%s\n',flat{:}); end
