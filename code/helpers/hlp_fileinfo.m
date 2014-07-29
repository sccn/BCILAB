function [version,inargs,outargs,h1,hash] = hlp_fileinfo(identifier,fname,whatinfo)
% Retuns (cached) information on a file (MD5 hash, input/output argument names).
% [Version,InArgs,OutArgs,H1Line,Hash] = hlp_fileinfo(Function,Filepath,WhatInfo)
%
% In:
%   Function : function handle or function name to look up
%   Filename : Optionally the path to the function (can disambiguate cases where Function is only a
%              name)
%   WhatInfo : Optionally the type of info to output as first output (can be any of the names under Out)
%
% Out:
%   Version : version string (if present) or hash (otherwise) for the function
%   InArgs  : names of input arguments
%   OutArgs : names of output arguments
%   H1Line  : the file's H1 line
%   Hash    : MD5 hash of the function
%
% Notes:
%   In deployed mode, this function relies on .m files being shipped separately with the binary,
%   using utl_whichfile.
%
% See also:
%   hlp_inoutargs
%
%				          		 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

persistent filedb;

if isdeployed
        
    % deployed mode must use a special mechanism (the m file source tree is maintained elsewhere)
    % no refreshing is necessary, as the function code is frozen
    if ~ischar(identifier)
        identifier = char(identifier); end
    try
        % look up the record for this file name
        fileentry = filedb.(identifier);
    catch %#ok<CTCH>
        % not found: look up file name (if you are not using BCILAB and get an error in the line
        % below this is because hlp_fileinfo will only work in deployed mode when a lot of
        % infrastructure is on the path)
        fname = utl_whichfile(identifier);
        fileentry.version = hlp_funcversion(identifier,fname);
        fileentry.md5 = hlp_cryptohash(fname,true);
        [fileentry.inargs,fileentry.outargs,fileentry.h1] = hlp_inoutargs(fname,identifier);
        filedb.(identifier) = fileentry;
    end
    inargs = fileentry.inargs;
    outargs = fileentry.outargs;
    version = fileentry.version;
    hash = fileentry.md5;
    h1 = fileentry.h1;    
    
else
    refresh_period = 1;  % refreshing our view of the file system at most every refresh_period seconds
    update = false;      % indicates whether the database entry for this file is valid or has to be (re)created
    
    % get the filename if necessary
    if nargin < 2 || isempty(fname)
        if ischar(identifier)
            fname = which(identifier);
        else
            funcinfo = functions(identifier);
            fname = funcinfo.file;
        end
    end
    
    if ~ischar(identifier)
        identifier = char(identifier); end
    
    try
        % look up the identifier from the DB
        identry = filedb.(identifier);
    catch %#ok<CTCH>
        % or create a blank identifier entry
        identry = struct('names',{{}},'entries',{{}});
        filedb.(identifier) = identry;
        update = true;
    end
    
    try
        % look up the file record for this identifier...
        idx = strcmp(identry.names,fname);
        fileentry = identry.entries{idx};
    catch %#ok<CTCH>
        % or create a blank file entry
        idx = length(identry.names)+1;
        identry.names{idx} = fname;
        identry.entries{idx} = [];
        fileentry = [];
        update = true;
    end
    
    % check if we need to recreate or update with data from the OS
    if update || (toc(fileentry.data.last_check) > refresh_period)
        filedata = dir(fname);
        % check if the file was actually changed
        if ~update && strcmp(fileentry.data.date,filedata.date)
            % remember when we last checked the OS record
            filedb.(identifier).entries{idx}.data.last_check = tic;
        else
            % the record has to be renewed
            try
                filedata(1).version = hlp_funcversion(identifier,fname);
            catch e
                hlp_handleerror(e);
            end
            filedata.md5 = hlp_cryptohash(fname,true);
            [filedata.inargs,filedata.outargs,filedata.h1] = hlp_inoutargs(fname);
            filedata.last_check = tic;
            fileentry.data = filedata;
            identry.entries{idx} = fileentry;
            filedb.(identifier) = identry;
        end
    end
    
    inargs = fileentry.data.inargs;
    outargs = fileentry.data.outargs;
    hash = fileentry.data.md5;
    version = fileentry.data.version;
    h1 = fileentry.data.h1;    
end

% optionally output the desired piece of info as first output (for convenience)
if nargin>=3
    switch lower(whatinfo)
        case 'hash'
            version = hash;
        case {'h1','h1line'}
            version = h1;
        case 'inargs'
            version = inargs;
        case 'outargs'
            version = outargs;
        case 'version'
            % nothing to do
        otherwise
            error(['Unsupported file info requested: ' whatinfo]);
    end
end
