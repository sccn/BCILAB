function io_save(fname, varargin)
% Like MATLAB's builtin save(), except that the fname can be platform-independent.
% io_save(Filename, Arguments...)
%
% Also, additional -options are supported, which offer additional fault-tolerance and automation,
% especially for long-running batch scripts.
%
% In:
%   Filename    :   Platform-Independent file name. If no ending is given, .mat will be assumed.
%                   If the ending is .sto, then BCILAB's Storage file format will be used (10-50x as
%                   fast as .mat).
%
%   Arguments...:   optional arguments to save; in addition to what is supported by save(), the following options are allowed:
%                   -makedirs   : make directories (recursively), if possible
%                   -retryinput : ask the user whether he/she wants to retry creating the directory or the file, if unsuccessful (using input())
%                   -attributes attribs : specify file/directory attributes; these are the second or second and third argument to be passed to 
%                                         the fileattrib() function, separated by a comma, to set attributes for the current user or specific user groups, 
%                                         e.g., '+w' or '+w','a'.
%                   -prunehandles : prune unreferenced variables from anonymous function handles; these are normally not used but saved anyway since cases
%                                   can be constructed in which they are actually needed by the function (e.g. using evalin('caller',...))
%                                   since these can be extremely large, they can easily cause save() to fail
%                   -serialized : compress the data using hlp_serialize
%                                 (when loading with io_load, the data will be automatically uncompressed,
%                                 not so when using load directly)
%
% Example:
%   % save the variables a,b, and c to a .mat file
%   io_save('store:/myoutput.mat', 'a','b','c');
%
%   % as before, but try to create directories if necessary
%   io_save('store:/myoutput.mat', 'a','b','c','-makedirs');
%
%   % like before, but this time prompt the user to re-try if an error occurs
%   io_save('store:/myoutput.mat', 'a','b','c','-makedirs','-retryinput');
%
%   % like before, but this time make the created file writable for all users
%   io_save('store:/myoutput.mat', 'a','b','c','-attributes','''+w'',''a''');
%
% See also:
%   save, io_load
%                     
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-21


% translate the file name
fname = env_translatepath(fname);
if length(fname)<4 || ~any(fname(end-3:end-1)=='.')
    fname = [fname '.mat']; end

% parse arguments...
saveargs = {'-v7.3'}; % args to save()
variables = {};       % names of variables to save (subset of saveargs)
makedirs = false;     % whether to create directories
retryinput = false;   % whether to ask the user for a retry
prunehandles = false; % whether to prune unreferenced workspace variables from function handles in arguments
serialized = false;   % whether to write variables as serialized
nooverwrite = false;  % whether to disable overwriting existing files
fileattriblist = [];  % the list of attributes for files & directories
i = 1;
while i <= length(varargin)
    arg = varargin{i};
    if ischar(arg) && ~isempty(arg) && arg(1) == '-'
        % we are dealing with an option argument
        switch(strtrim(arg))
            case '-makedirs'
                makedirs = true;
            case '-retryinput'
                retryinput = true;
            case '-attributes'
                fileattriblist = varargin{i+1};
                i = i+1; % skip next arg, too
            case '-prunehandles'
                prunehandles = true;
            case {'-serialized','-serialize'}
                serialized = true;
            case '-nooverwrite'
                nooverwrite = true;
            otherwise
                % regular argument; append it
                saveargs{end+1} = arg;
        end
    else
        % append the argument
        saveargs{end+1} = arg;
        variables{end+1} = arg;
    end
    i = i+1; % go to next argument
end

% move variables into a local workspace for optional processing
if isempty(variables)
    % no variables listed: obtain all
    varinfos = evalin('caller','whos');
    variables = {varinfos.name};
    saveargs = [saveargs variables];
end
for k=1:length(variables)
    wkspace.(variables{k}) = evalin('caller',variables{k}); end

% prune their handles
if prunehandles
    for k=1:length(variables)
        wkspace.(variables{k}) = utl_prune_handles(wkspace.(variables{k})); end
end

% serialize the variables
if serialized
    for k=1:length(variables)
        wkspace.(variables{k}) = hlp_serialize(wkspace.(variables{k})); end
    wkspace.is_serialized__ = ['variables were serialized using ' hlp_funcversion(@hlp_serialize) '; use hlp_deserialize (see File Exchange) to return them to their original form.'];
    saveargs{end+1} = 'is_serialized__';
end

% reformat the save-args into a string
tmp = [saveargs' repmat({' '},length(saveargs),1)]'; 
saveargs = [tmp{:}];

% reformat the fileattriblist into a cell array...
fileattriblist = eval(['{' fileattriblist '}']);

while 1    
    if makedirs
        % first ensure that the directory exists
        while 1
            try
                io_mkdirs(fname,fileattriblist);
                break;
            catch
                warning('BCILAB:io_save:cannot_create_directory',['Error creating target directory ' fname '. Please review your permissions.']);
                if retryinput
                    if strcmp(input('Would you like to retry? [yes/no]: ','s'),'no')
                        disp('Did not save the data.');
                        break;
                    end
                else
                    break;
                end
            end
        end
    end
    
    % we have copied the variables to be saved into the 'wkspace' variable; now save
    try
        if nooverwrite && exist(fname,'file')
            return; end
        if length(fname)>4 && strcmpi(fname(end-3:end),'.sto')
            % save as .sto file (fast)
            bytes = hlp_serialize(wkspace);
            f = fopen(fname,'w');
            fwrite(f,bytes);
            fclose(f);
        else            
            % save as .mat file (slow)
            save_workspace(fname,saveargs,wkspace);
        end
        if ~isempty(fileattriblist)
            % try to set file fileattriblist...
            if ~isempty(fileattriblist)
                warning off MATLAB:FILEATTRIB:SyntaxWarning
                fileattrib(fname,fileattriblist{:}); 
            end            
        end
        break;
    catch e
        warning('BCILAB:io_save:cannot_save_file',['Error saving target file ' fname '. Please review your permissions: ' e.message]);
        if retryinput
            if strcmp(input('Would you like to retry? [yes/no]','s'),'no') 
                disp('Did not save the data.');
                break;
            end
        else
            break;
        end
    end
end

% save, given a workspace of variables
function save_workspace(fname__,saveargs__,workspace__)
for fn__=fieldnames(workspace__)'
    eval([fn__{1} ' = workspace__.(fn__{1});']); end
clear workspace__ fn__;
eval(['save ''' fname__ ''' ' saveargs__ ';']);
