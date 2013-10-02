function collected = hlp_collect_datasets(directory,varargin)
% Find all loadable data sets in some directory tree and collect arbitary properties.
% Collected = hlp_collect_datasets(Directory, Options...)
%
% In:
%   Directory    : the directory to be searched for eeg data sets
%
%   Options...   : optional name-value pairs; possible names are:
%                   'pattern': regular expression pattern of file paths considered (default: regexptranslate('wildcard','*.set'))
%                   'nopattern': cell array of regular expression patterns (for files/directories) to exclude from scan (default: {})
%                   'checkset': whether to run an eeg_checkset on the data (default: 1)
%                   'nowarnings': exclude files that give a warning when loading or running eeg_checkset (default: 1)
%                   'nodialogs': exclude files that create dialogs when loading or running eeg_checkset (default: 1)
%                                note: this requires a certain dependency directory (and uses env_translatepath to find it)
%                   'conditions': cell array of functions that check custom conditions on the data (e.g., @(EEG) ~isempty(EEG.icawinv)) (default: [])
%                   'maxsize': maximum filesize considered, in bytes (default: 5*2^20 = 5MB)
%                   'maxtime': maximum allowed file processing time, in seconds (default: Inf)
%                   'maxnumber': maximum number of entries returned (default: Inf)
%                   'collect': a function of the EEG set (and optionally file path); the function's return values are collected
%                              for every admissible data set (default: @(EEG,path) path)
%                   'fileconditions': cell array of functions that check custom conditions on the filename (e.g. @(path,name,ext) exist([path filesep name '.xxx'],'file'))
%                                     (default: @(path,name,ext) exist([path filesep name '.fdt'],'file') || exist([path filesep name '.dat'],'file'))
%
% Out:
%   Collected : cell array of admissible file paths (or desired contents, if 'collect' was specified)
%
% Notes:
%   If the function is terminated prematurely, the global variables collected_so_far and num_collected_so_var give the current data:
%   global collected_so_far num_collected_so_far; mydata = collected_so_far(1:num_collected_so_far)
%
% Examples:
%   % collect .icawinv and .chanlocs fields of all data sets for which more than half of the channels have locations
%   collected = hlp_collect_datasets('/data/projects', ...
%       'conditions',@(EEG) ~isempty(EEG.icawinv) && mean(cellfun('isempty',{EEG.chanlocs.X})) < 0.5, ...
%       'collect',@(EEG,path){path,EEG.icawinv,EEG.chanlocs,EEG.icachansind})
%
%   % as before, but consider files that are up to 100MB large
%   collected = hlp_collect_datasets('/data/projects', 'maxsize',100*2^20, ...
%       'conditions',@(EEG) ~isempty(EEG.icawinv) && mean(cellfun('isempty',{EEG.chanlocs.X})) < 0.5, ...
%       'collect',@(EEG,path){path,EEG.icawinv,EEG.chanlocs,EEG.icachansind})
%
%   % like before, but this time do not scan directories which contain the strings 'christian' and 'duann'
%   collected = hlp_collect_datasets('/data/projects', 'nopattern',{'christian','duann'}, ...
%       'conditions',@(EEG) ~isempty(EEG.icawinv) && mean(cellfun('isempty',{EEG.chanlocs.X})) < 0.5, ...
%       'collect',@(EEG,path){path,EEG.icawinv,EEG.chanlocs,EEG.icachansind})
%
%   % like before, but stop after having gathered 200 data sets
%   collected = hlp_collect_datasets('/data/projects', 'maxnumber',200, ...
%       'conditions',@(EEG) ~isempty(EEG.icawinv) && mean(cellfun('isempty',{EEG.chanlocs.X})) < 0.5, ...
%       'collect',@(EEG,path){path,EEG.icawinv,EEG.chanlocs,EEG.icachansind})
%
%                               Laura Froelich, Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-10

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

global collected_so_far num_collected_so_far;                                % collected_so_far(1:numcollected_so_far) is what we have so far...

if length(varargin)==1 && isstruct(varargin{1})
    opts = varargin{1};                                                     % opts is already a struct: fast path
else
    opts = hlp_varargin2struct(varargin, 'pattern',regexptranslate('wildcard','*.set'), 'nopattern',{}, 'checkset',1, 'nowarnings',1 ,'nodialogs',1, ...
        'maxsize',5*2^20, 'maxtime',Inf, 'maxnumber',Inf, 'conditions',{}, 'collect',@(EEG,path) path, ....
        'fileconditions',@(path,name,ext) exist([path filesep  name '.fdt'],'file') || exist([path filesep name '.dat'],'file'));
    if ~iscell(opts.conditions)
        opts.conditions = {opts.conditions}; end
    if ~iscell(opts.nopattern)
        opts.nopattern = {opts.nopattern}; end
    if ~iscell(opts.fileconditions)
        opts.fileconditions = {opts.fileconditions}; end
    if opts.nodialogs
        % the nodialogs option unfortunately requires some MATLAB hacking...
        if isdeployed
            warning('Dialogs cannot be disabled in deployed mode.'); end %#ok<WNTAG>
        if hlp_matlab_version < 706
            warning('BCILAB:hlp_collect_datasets:nodialogs_pathmess','With MATLABs older than 2008a, using the ''nodialogs'' option will leave you in another directory after running hlp_collect_datasets.'); 
        end
        % need to be able to call this recursively after we cd'd to another path
        thisdir = fileparts(mfilename('fullpath'));
        addpath(thisdir);
        olddir = pwd;
        if exist('env_translatepath','file')
            % if this function is present, we can directly resolve the location of the folder with the disabled dialogs
            cd(env_translatepath('dependencies:/disabled_dialogs'));
        elseif exist([thisdir filesep 'private' filesep 'disabled_dialogs'],'dir')
            cd([thisdir filesep 'private' filesep 'disabled_dialogs']);
        else
            warning('BCILAB:hlp_collect_datasets:nodialogs_dysfunctional','Please create a directory called ''private'' in %s and put the "disabled_dialogs" folder in there (which you should find in your BCILAB dependencies directory).',thisdir);
        end
        go_back = onCleanup(@()cd(olddir));
    end
    collected_so_far = {};
    num_collected_so_far = 0;
end

disp(['entering ' directory '...']);
collected = [];
topfiles = dir(directory);
topfiles = topfiles([topfiles.bytes] <= opts.maxsize);                      % discard too large files
for it = {topfiles.name}                                                    % for each admissible dir entry
    item = it{1};
    whole_path = [directory filesep item];
    if any(cellfun(@(x) ~isempty(regexp(whole_path,x,'once')), opts.nopattern)) % discard disallowed patterns
        disp(['skipping ' whole_path '...']);
        continue;
    end
    if isdir(whole_path)
        if ~isempty(item) && item(1)~='.'                                   % discard self & parent paths, and hidden paths
            collected = [collected hlp_collect_datasets(whole_path, opts)]; end %#ok<AGROW> % recurse...
    elseif regexp(whole_path,opts.pattern)                                  % discard non-matching files
        try
            fprintf(['testing ' whole_path '... ']);
            [path,name,ext] = fileparts(whole_path);
            for cond = opts.fileconditions                                  % check for file name conditions
                if cond{1}(path,name,ext)
                    % succeeded
                else
                    error('file name condition violated');
                end
            end
            t0 = tic;                                                        % measure processing time
            [conout,data] = evalc(sprintf('pop_loadset(''filename'',''%s'', ''filepath'',''%s'', ''loadmode'',''info'', ''check'',''off'')',item,directory));
            if opts.nowarnings && ~isempty(strfind(lower(conout),'warning')) % discard files with load warnings
                error('loadset warning'); end
            if opts.checkset
                [conout,data] = evalc('eeg_checkset(data)');
                if opts.nowarnings && ~isempty(strfind(lower(conout),'warning')) % discard files with checkset warnings
                    error('checkset warning'); end
            end
            if toc(t0) >= opts.maxtime                                       % discard files that take too long to process
                error('processing time exceeded'); end
            for cond = opts.conditions                                       % check for additional conditions
                if cond{1}(data)
                    % succeeded
                else
                    error('data set condition violated');
                end
            end
            if nargin(opts.collect) == 2                                    % collect properties
                selection = opts.collect(data,whole_path);
            else
                selection = opts.collect(data);
            end
            collected{end+1} = selection; %#ok<AGROW>
            num_collected_so_far = num_collected_so_far+1;
            if length(collected_so_far) < num_collected_so_far
                collected_so_far{1+2*end} = []; end                         % grow this array in ~ constant time
            collected_so_far{num_collected_so_far} = selection;             % track results globally, too
            disp('included.');
        catch e
            disp(['excluded: ' e.message]);
        end
        if length(collected) >= opts.maxnumber
            break; end
    end
end
