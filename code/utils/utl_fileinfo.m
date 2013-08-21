function [version,inargs,outargs,h1,hash] = utl_fileinfo(fname,identifier)
% Retuns (cached) information on a file (MD5 hash, input/output argument names).
% [Version,InArgs,OutArgs,H1Line,Hash] = utl_fileinfo(Filename,Identifier)
%
% In:
%   Filename   : full pathname of an m-file
%   Identifier : MATLAB identifier of the file
%
% Out:
%   Version : version string (if present) or hash (otherwise) for the function
%   InArgs  : names of input arguments
%   OutArgs : names of output arguments
%   H1Line  : the file's H1 line
%   Hash    : MD5 hash of the function
%
% Notes:
%   In deployed mode, this function relies on .m files being shipped separately with the binary.
%   If the global cell string tracking.paths.toolboxes exists, its contents will be used as root 
%   source directories; otherwise, env_translatepath('functions:/') will be used.
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
    try
        % look up the record for this file name
        fileentry = filedb.(identifier);
    catch
        % not found: look up file name
        fname = utl_whichfile(identifier);
        fileentry.version = utl_funcversion(identifier);
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
    
    try
        % look up the identifier from the DB
        identry = filedb.(identifier);
    catch
        % or create a blank identifier entry
        identry = struct('names',{{}},'entries',{{}});
        filedb.(identifier) = identry;
        update = true;
    end
    
    try
        % look up the file record for this identifier...
        idx = strcmp(identry.names,fname);
        fileentry = identry.entries{idx};
    catch
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
            filedata.version = utl_funcversion(identifier);
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
