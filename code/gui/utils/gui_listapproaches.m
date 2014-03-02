function [list,flat,indexable] = gui_listapproaches()
% Compile a list of all approaches that are currently accessible
% [List,Flat,Indexable] = gui_listapproaches()
%
% The list is compiled from the workspace variables, the original paradigms (in the paradigms folder), and the files in the 'resources/approaches' directory...
%
% Out:
%   List : cell array of {'Groupname', {approach1,approach2,approach3, ...}, 'Groupname', {approach1,approach2,approach3}, ...}
%
%   Flat : pretty-printed cell-string array of List contents
%
%   Indexable : a flat cell array of approaches that can be indexed in accordance with Flat
%
% Notes:
%   if no output is taken, the function displays the Flat list.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

global tracking;

list = {};

% --- aggregate all paradigms ---

paradigm_locations = { ...
    {'Built-in Paradigms', {'functions:/paradigms'}}, ...
    {'User Paradigms', {'home:/.bcilab/code/paradigms'}}, ...
    {'Experimental Paradigms',quickif(tracking.gui.show_experimental,{'functions:/paradigms/in_development'},{})}};
if ~isempty(tracking.paths.private_path)
    paradigm_locations = [paradigm_locations { ...
        {'Private Paradigms', {'private:/code/paradigms'}}, ...
        {'Private Experimental Paradigms',quickif(tracking.gui.show_experimental,{'private:/code/paradigms/in_development'},{})}}];
end

for loc=1:length(paradigm_locations)
    if ~isempty(paradigm_locations{loc}{2})
        para_paths = {};
        for p = paradigm_locations{loc}{2}
            rootpath = env_translatepath(p{1});
            files = dir([rootpath filesep 'Paradigm*.m']);
            para_paths = [para_paths cellfun(@(n)[rootpath filesep n],{files.name},'UniformOutput',false)];
        end
        if ~isempty(para_paths)
            list = [list {paradigm_locations{loc}{1}, list_paradigms_in_paths(para_paths)}]; end
    end    
end

% -- get list of paradigm names ---
paradigm_names = {};
for l=2:2:length(list)
    for p = list{l}
        paradigm_names{end+1} = p{1}.paradigm(9:end); end
end

% --- aggregate all approaches in the workspace ---
approaches = {};

vars = evalin('base','whos');
for v=find(strcmp({vars.class},'struct'))
    var = evalin('base',vars(v).name);
    if isfield(var,{'paradigm','parameters'})
        % found one: list it
        if isfield(var,'name')
            var.name = [var.name ' (' vars(v).name ')']; 
        else
            var.name = vars(v).name; 
        end
        approaches{end+1} = var;
    end
end

for v=find(strcmp({vars.class},'cell'))
    var = evalin('base',vars(v).name);
    if length(var)>1 && ischar(var{1}) && ((strncmp(var{1},'Paradigm',8) && length(var{1})>8 && any(strcmp(var{1}(9:end),paradigm_names))) || any(strcmp(var{1},paradigm_names)))
        approaches{end+1} = struct('paradigm',quickif(strncmp(var{1},'Paradigm',8),var{1},['Paradigm',var{1}]),'name',vars(v).name,'parameters',{var(2:end)},'description',''); end
end

list = [list {'From Workspace', approaches}];


% --- aggregate all approaches in the approaches directory (and the user's home directory, too) ---
approaches = {};
approach_dirs = {env_translatepath('home:/.bcilab/approaches'),env_translatepath('resources:/approaches')};
if ~isempty(tracking.paths.private_path)
    approach_dirs = [approach_dirs {env_translatepath('private:/approaches')}]; end
for d = approach_dirs
    approach_dir = d{1};
    files = dir([approach_dir filesep '*.apr']);
    filenames = {files.name};

    for f=1:length(filenames)
        % load the file
        try
            contents = io_load([approach_dir filesep filenames{f}],'-mat');
            for fn=fieldnames(contents)'
                var = contents.(fn{1});
                if isfield(var,{'paradigm','parameters'})
                    % found one: list it
                    if isfield(var,'name')
                        var.name = [var.name ' (' fn{1} ')'];
                    else
                        var.name = fn{1};
                    end
                    approaches{end+1} = var;
                end
            end
        catch e
            disp(['Could not load approach ' filenames{f} ' from disk.']);
        end
    end
end

list = [list {'From Disk', approaches}];

[flat,indexable] = gui_print_grouplist(list,@(x)x.name);

% print it
if nargout == 0
    fprintf('%s\n',flat{:}); end



function approaches = list_paradigms_in_paths(para_paths)
persistent cached_paradigms;
approaches = {};
for p=1:length(para_paths)
    [paradir,filename,ext] = fileparts(para_paths{p});
    filename = [filename ext];
    tag = filename(1:end-2);
    try
        % try lookup from a tiny cache of approach descriptions        
        appr = cached_paradigms.(tag);
    catch
        % otherwise generate it
        try
            % read get the help text
            % (we don't even try to use help() as it doesn't work in deployed mode)
            f = fopen([paradir filesep filename],'r');
            help_lines = {};
            while 1
                line = fgetl(f);
                % end of file?
                if ~ischar(line)
                    break; end
                line = strtrim(line);
                % function declaration line?
                if strncmp(line,'classdef',8)
                    continue; end
                % non-comment line?
                if isempty(line) || line(1) ~= '%'
                    break; end
                % methods line?
                if strncmp(line,'methods',8)
                    break; end
                % strip leading comment chars
                while ~isempty(line) && line(1) == '%'
                    line(1) = []; end
                % append
                help_lines{end+1} = [' ' line];
            end
            fclose(f);
        catch
            try fclose(f); catch,end
            error('Cannot read help text for file %s.',filename);
        end
        if isempty(help_lines)
            help_lines = {'help for file: ',[paradir filesep filename]}; end
        
        % extract the description
        name_line = find(~cellfun('isempty',strfind(help_lines, 'Name:')));
        desc_lines = help_lines(1:name_line-2);
        % deblank the lines
        for l=1:length(desc_lines)
            desc_lines{l} = deblank(desc_lines{l}); end
        % find out if we can dedent the entire help...
        min_indent = Inf;
        for l=1:length(desc_lines)
            dl = desc_lines{l};
            if ~isempty(dl)
                min_indent = min(min_indent,length(dl) - length(strtrim(dl))); end
        end
        if min_indent > 0
            for l=1:length(desc_lines)
                if ~isempty(desc_lines{l})
                    desc_lines{l} = desc_lines{l}(min_indent+1:end); end
            end
        end
        % reformat line breaks and spaces
        for l=1:length(desc_lines)
            dl = desc_lines{l};
            if isempty(dl)
                dl = [dl sprintf('\n\n')];
            elseif length(strtrim(dl)) < length(dl) || strcmp(dl,'References:')
                dl = [dl sprintf('\n')];
            else
                dl = [dl ' '];
            end
            desc_lines{l} = dl;
        end
        desc = [desc_lines{:}];
        if ~isempty(desc)
            desc = strrep(desc,'  ',' '); end
        % extract the name
        if isempty(name_line)
            name = filename(1:end-2);
        else
            name = [strtrim(help_lines{name_line+1}) ' (' filename(1:end-2) ')'];
        end
        appr = struct('paradigm',filename(1:end-2), 'parameters',{{}}, 'description',desc, 'name',name);
        % and store it for the next time
        cached_paradigms.(tag) = appr;
    end
    if isempty(strfind(appr.name,'(abstract'))
        approaches{end+1} = appr; end
end
