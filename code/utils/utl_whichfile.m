function filename = utl_whichfile(identifier)
% Find the .m file which defines the given function identifier.
% Filename = utl_whichfile(Identifier)
%
% This is a special-purpose alternative to the builtin 'which', for use in deployed situations.
%
% In these cases, 'which' returns the filename of only an encrypted version of the defining
% function, so it cannot be used to access the m code. This function assumes that the source code
% for each function ships with the release, and returns the file name for a given identifier.
%
% In:
%   Identifier : identifier of some function
%
% Out:
%   Filename : .m file that contains the code which implements the
%              identified function
%
% Notes:
%   This function does not have all the facilities of which(), and assumes that there is one unique
%   source file for every looked up identifier.
%
%   If the global cell string tracking.paths.toolboxes exists, its contents will be used as root
%   directories; otherwise, env_translatepath('functions:/') will be used.
%
% See also:
%   which
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-04

persistent filenames; % table of m file names, indexed by identifier

try
    % look up the record for this file name
    filename = filenames.(identifier);
catch
    if isempty(filenames)
        % file table is empty (first-time run): build it recursively from root directories
        global tracking; %#ok<TLEV>
        if isfield(tracking,'paths') && isfield(tracking.paths,'toolboxes')
            root_dirs = tracking.paths.toolboxes;
        else
            root_dirs = {env_translatepath('functions:/')};
        end
        for d=1:length(root_dirs)
            populate_records(root_dirs{d}); end
    else
        error('No record for the file %s was found. Make sure it is in the searched directories.',identifier);
    end
    try
        filename = filenames.(identifier);
    catch
        error('No record for the file %s was found. Make sure it is in the searched directories.',identifier);
    end
end

    % recursively populate file table from an initial directory.
    function populate_records(path)
        if path(end) == filesep
            path(end) = []; end
        for info = dir(path)'
            if info.isdir
                % found a directory: recurse
                if ~any(strcmp(info.name,{'.','..'}))
                    populate_records([path filesep info.name]); end
            else
                % found a file --? add it to the table if an .m file
                if length(info.name) > 2 && strcmp(info.name(end-1:end),'.m')                    
                    filenames.(info.name(1:end-2)) = [path filesep info.name]; end
            end
        end
    end

end