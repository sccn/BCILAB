function run_script(script,runinbase)
% Run a script in the caller's workspace.
%
% In:
%   Script: name of the script
%
%   BaseWorkspace : if true, run the script in the 'base' workspace (default: false)
%
% Notes:
%   Certain functions that refer to the name of the currently executing file may yield unexpected
%   results when using this function (so far except for mfilename('fullpath'), which is supported).
%
% See also:
%   run
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-23

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

if exist('runinbase','var') && runinbase
    runworkspace = 'base';
else
    runworkspace = 'caller';
end


% correct platform-dependent paths
script = strrep(strrep(script,'\',filesep),'/',filesep);

% expand relative paths
if ~any(script == filesep)
    script = [pwd filesep script]; end

% add default extension
if exist([script '.m'],'file')
    script = [script '.m']; end

% check for existence
if ~exist(script,'file')
    error('run_script:file_not_found','Can''t find script %.',script); end

% try to read line by line
try
    lines = {};
    f = fopen(script,'r');
    while 1
        l = fgetl(f);
        if ~ischar(l)
            break; end
        lines{end+1} = [l 10];
    end
    fclose(f);
catch e
    try fclose(f); catch,end
    error('run_script:cannot_read_file','Can''t read script %s; Error message: %s',script,e.message);
end

% turn into one string that can be passed to evalin
% this requires that newlines and ellipses are properly handled...
evalstr = char([lines{:}]);
comment_flag = false;
string_flag = false;
ellipsis_flag = false;
spaceout = false(1,length(evalstr));   % this mask indicates where we can substitute irreversibly by whitespace characters...
for k=1:length(evalstr)
    if ellipsis_flag
        % everything that follows an ellipsis will be spaced out (including the subsequent newline that resets it)
        spaceout(k) = true; end    
    switch evalstr(k)
        case '''' % quotes
            % flip str flag, unless in comment
            if ~comment_flag
                string_flag = ~string_flag; end
        case 10 % newline
            % reset bracket level, unless in ellipsis
            % reset comment flag, str flag and ellipsis flag
            comment_flag = false;
            string_flag = false;
            ellipsis_flag = false;
        case '%' % comment character
            % if not in str, switch on comment flag
            if ~string_flag
                comment_flag = true; end
        case '.' % potential ellipsis character
            % if not in comment nor in str, turn on ellipsis and comment
            if ~string_flag && ~comment_flag && k>2 && strcmp(evalstr(k-2:k),'...')
                ellipsis_flag = true;
                comment_flag = true;
                % we want to replace the ellipsis and everything that follows up to and including the next newline
                spaceout(k-2:k) = true;
            end
    end
end


% space out the ellipses
evalstr(spaceout) = ' ';

% replace mfilename('fullpath') with '/path/to/script'
evalstr = strrep(evalstr,'mfilename(''fullpath'')',['''' script '''']);

% cd into the script directory, but make sure that we go back once completed
olddir = pwd;
cd(fileparts(script));

if exist('onCleanup','file')
    % use onCleanup to get back to the old directory
    go_back = onCleanup(@()cd(olddir));
    % evaluate the script
    evalin(runworkspace,evalstr);
else
    % manually go back to the old directory, even in case of an exception
    % note: does not work if Ctrl+C is pressed during script execution
    try
        % evaluate in the caller's workspace
        evalin(runworkspace,evalstr);
        cd(olddir);
    catch e
        cd(olddir);
        rethrow(e);
    end
end