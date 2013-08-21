function result = hlp_config(configname, operation, varargin)
% helper function to process human-readable config scripts.
% Result = hlp_config(FileName,Operation,VariableName,Value,NVPs...)
%
% Config scripts consist of assignments of the form name = value; to set configuration options. In
% addition, there may be any type of comments, conditional control flow, etc - e.g., setting certain 
% values on some platforms and others on others. This function allows to get or set the value 
% assigned to a variable in the place of the script where it is actually assigned on the current 
% platform. Note that the respective variable has to be already in the config file for this function 
% to work.
%
% In:
%   FileName : name of the configuration file to process
%
%   Operation : operation to perform on the config file
%               'get' : get the currently defined value of a given variable
%               'set' : replace the current defintion of a given variable
%
%   VariableName : name of the variable to be affected (must be a MATLAB identifier)
%
%   Value : the new value to be assigned, if the operation is 'set', as a string
%           note that most data structures can be converted into a string via hlp_tostring
%
%   NVPs... : list of further name-value pairs, where each name denotes a config variables and the subsequent
%             value is the string expression that should be written into the config file. It is 
%             generally a good idea to use hlp_tostring() to turn a data structure into such a string
%             representation.
%
% Out:
%   Result : the current value of the variable of interest, when using the 'get'
%            operation
%
% Notes:
%   There can be multiple successive variable name / value pairs for the set mode.
%   If an error occurs during a set operation, any changes will be rolled back.
% 
% Examples:
%   % read out the value of the 'data' config variable from a config file
%   data = hlp_config('/home/christian/myconfig.m','get','data')
%
%   % override the values of the 'files' and 'capacity' config variables in the given config script
%   hlp_config('/home/christian/myconfig.m', 'set', 'files',myfiles, 'capacity',1000)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-19

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

if ~exist(configname,'file')
    error('hlp_config:file_not_found','The specified config file was not found.'); end

switch operation
    case 'get'
        varname = varargin{1};
        if ~isvarname(varname)
            error('hlp_config:bad_varname','The variable name must be a valid MATLAB identifier.'); end
        % get the currently defined value of a variable...
        result = get_value(configname,varname);
    case 'set'
        backupfile = [];
        try
            % apply first assignment
            backupfile = set_value(configname,varargin{1},varargin{2},true);
            for k = 4:2:length(varargin)
                % apply all other assignments
                set_value(configname,varargin{k-1},varargin{k},false); end
        catch e
            % got an error; roll back changes if necessary
            if ~isempty(backupfile)
                try
                    movefile(backupfile,configname); 
                catch
                    disp(['Could not roll back changes. You can manually revert changes by replacing ' configname ' by ' backupfile '.']);
                end
            end
            rethrow(e);
        end
    otherwise
        error('hlp_config:unsupported_option','Unsupported config operation.');
end


% run the given config script and obtain the current value of the given variable...
function res = get_value(filename__,varname__)
try
    run_script(filename__);
catch e
    error('hlp_config:erroneous_file',['The config file is erroneous; Error message: ' e.message]);
end
if ~exist(varname__,'var')
    error('hlp_config:var_not_found','The variable is not being defined in the config file.'); end
res = eval(varname__);


function backup_name = set_value(filename,varname,newvalue,makebackup)
backup_name = [];
if ~exist(filename,'file')
    error('hlp_config:file_not_found','The config file was not found.'); end
if ~isvarname(varname)
    error('hlp_config:incorrect_value','The variable name must be a valid MATLAB identifier.'); end
if ~ischar(newvalue)
    error('hlp_config:incorrect_value','The value to be assigned must be given as a string.'); end
try
    % read the config file contents
    contents = {};
    f = fopen(filename,'r');
    while 1
        l = fgetl(f);
        if ~ischar(l)
            break; end
        contents{end+1} = [l 10];
    end
    fclose(f);
    % turn it into one str
    contents = [contents{:}];
catch e
    try fclose(f); catch,end
    error('hlp_config:cannot_read_config',['Cannot read the config file; Error message: ' e.message]);
end

% now check if the file is actually writable
try
    f = fopen(filename,'r+');
    if f ~= -1
        fclose(f);
    else
        error('hlp_config:permissions_error','Could not update the config file %s. Please check file permissions and try again.',filename);
    end
catch
    error('hlp_config:permissions_error','Could not update the config file %s. Please check file permissions and try again.',filename);
end

% temporarily replace stray semicolons by a special character and contract ellipses,
% so that the subsequent assignment regex matching will not get derailed)
evalstr = contents;
comment_flag = false;
string_flag = false;
bracket_level = 0;
ellipsis_flag = false;
substitute = false(1,length(evalstr)); % this mask indicates where we have to subsitute reversibly by special characters
spaceout = false(1,length(evalstr));   % this mask indicates where we can substitute irreversibly by whitespace characters...
for k=1:length(evalstr)
    if ellipsis_flag
        % everything that follows an ellipsis will be spaced out (including the subsequent newline that resets it)
        spaceout(k) = true; end    
    switch evalstr(k)
        case ';' % semicolon
            % in strs, brackets or comments: indicate need for substitution
            if string_flag || bracket_level>0 || comment_flag
                substitute(k) = true; end
        case '''' % quotes
            % flip str flag, unless in comment
            if ~comment_flag
                string_flag = ~string_flag; end
        case 10 % newline
            % reset bracket level, unless in ellipsis
            if ~ellipsis_flag
                bracket_level = 0; end
            % reset comment flag, str flag and ellipsis flag
            comment_flag = false;
            string_flag = false;
            ellipsis_flag = false;
        case {'[','{'} % opening array bracket
            % if not in str nor comment, increase bracket level
            if ~string_flag && ~comment_flag
                bracket_level = bracket_level+1; end
        case {']','}'} % closing array bracket
            % if not in str nor comment, decrease bracket level
            if ~string_flag && ~comment_flag
                bracket_level = bracket_level-1; end
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
% replace the characters that need to be substituted (by the bell character)
evalstr(substitute) = 7;
evalstr(spaceout) = ' ';
% replace all assignments of the form "varname = *;" by "varname{end+1} = num;"
[starts,ends] = regexp(evalstr,[varname '\s*=[^;\n]*;']);
for k=length(starts):-1:1
    evalstr = [evalstr(1:starts(k)-1) varname '{end+1} = struct(''assignment'',' num2str(k) ');' evalstr(ends(k)+1:end)]; end
% add initial assignment
evalstr = [sprintf('%s = {};\n',varname) evalstr];
% back-substitute the special character by semicolons
evalstr(evalstr==7) = ';';

% evaluate contents and get the matching assignment id's
ids = run_protected(evalstr,varname);

% check validity of the updated value, and of the updated config file
try
    % check if the value str can in fact be evaluated
    newvalue_eval = eval(newvalue);
catch
    error('hlp_config:incorrect_value','The value "%s" (to be assigned to variable "%s") cannot be evaluated properly. Note that, for example, string values need to be quoted.',newvalue,varname);
end
% evaluate the original config script and record the full variable assignment
[dummy,wspace_old] = run_protected(contents); %#ok<ASGLU>
% splice the new value into the config file contents, for the last assignment in ids
id = ids{end}.assignment;
contents = [contents(1:starts(id)-1) varname ' = ' newvalue ';' contents(ends(id)+1:end)];
% evaluate the new config script and record the full variable assignment
[dummy,wspace_new] = run_protected(contents); %#ok<ASGLU>
% make sure that the only thing that has changed is the assignment to the variable of interest
wspace_old.(varname) = newvalue_eval;
if ~isequalwithequalnans(wspace_old,wspace_new)
    error('hlp_config:update_failed','The config file can not be properly updated.'); end

% apparently, everything went well, except for the following possibilities
%  * the newly assigned value makes no sense (--> usage error)
%  * the settings were changed for unanticipated platforms (--> this needs to be documented properly)
if makebackup
    try
        % make a backup of the original config file using a fresh name (.bak00X)
        [p,n,x] = fileparts(filename);
        files = dir([p filesep n '*.bak*']);
        backup_numbers = cellfun(@(n)str2num(n(end-2:end)),{files.name},'UniformOutput',false);
        backup_numbers = [backup_numbers{:}];
        if ~isempty(backup_numbers)
            new_number = 1 + max(backup_numbers);
        else
            new_number = 1;
        end
        backup_name = [p filesep n '.bak' sprintf('%03i',new_number)];
        copyfile(filename,backup_name);
        % set read permissions
        warning off MATLAB:FILEATTRIB:SyntaxWarning
        fileattrib(backup_name,'+w','a');
    catch
        error('hlp_config:permissions_error','Could not create a backup of the original config file %s. Please check file permissions and try again.',filename);
    end
end
    
% split the contents into lines again
contents = strsplit(contents,10);
try
    % re-create the file, line by line
    f = fopen(filename,'w+');
    for k=1:length(contents)
        fwrite(f,contents{k});
        fprintf(f,'\n');
    end
    fclose(f);
    % set file attributes
    warning off MATLAB:FILEATTRIB:SyntaxWarning
    fileattrib(filename,'+w','a');
catch
    try fclose(f); catch,end
    error('hlp_config:permissions_error','Could not override the config file %s. Please check file permissions and try again.',filename);
end



% run the given config script and obtain the current value of the given variable...
function [res,wspace] = run_protected(code__,varname__)
try
    eval(code__);
    % collect all variables into a workspace struct
    infos = whos();
    for n = {infos.name}
        if ~any(strcmp(n{1},{'code__','varname__'}))
            wspace.(n{1}) = eval(n{1}); end
    end
    if exist('varname__','var')
        % if a specific variable was to be inspected...
        res = eval(varname__);
        if ~iscell(res) || length(res) < 1 || ~all(cellfun('isclass',res,'struct')) || ~all(cellfun(@(x)isfield(x,'assignment'),res))
            error('Not all assignments to the variable were correctly identified.'); end
    else
        res = [];
    end
catch e
    error('hlp_config:update_error',['The config file could not be parsed (probably it is ill-formed); Debug message: ' e.message]);
end



% split a string without fusing delimiters (unlike hlp_split)
function strs = strsplit(str, delim)
idx = strfind(str, delim);
strs = cell(numel(idx)+1, 1);
idx = [0 idx numel(str)+1];
for k = 2:numel(idx)
    strs{k-1} = str(idx(k-1)+1:idx(k)-1); end



% for old MATLABs that can't properly move files...
function movefile(src,dst)
try
    builtin('movefile',src,dst);    
catch e
    if any([src dst]=='$') && hlp_matlab_version <= 705
        if ispc
            [errcode,text] = system(sprintf('move ''%s'' ''%s''',src,dst)); %#ok<NASGU>
        else
            [errcode,text] = system(sprintf('mv ''%s'' ''%s''',src,dst)); %#ok<NASGU>
        end
        if errcode
            error('Failed to move %s to %s.',src,dst); end
    else
        rethrow(e);
    end
end