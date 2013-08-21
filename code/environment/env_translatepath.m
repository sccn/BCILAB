function filename = env_translatepath(filename)
% Translates platform-independent directories into a system-specific directories.
% SystemPath = env_translatepath(IndependentPath)
%
% BCILAB supports platform-independent paths for all its scripts and IO functions, which allows for
% script portability from system to system. It is especially important when data paths are mounted
% at different locations, depending on access mode and operating system. A side effect of portable
% path names is that there is one unique expression which computes any given data set, such as
% "flt_iir(flt_reref((flt_resample(io_loadset('data:/Projects/Test/test1.vhdr'),[],200))),[],[4 6 25
% 30])", and this in turn allows to share the same data set caches (which are indexed by expression)
% across machines and operating systems, minimizing redundant computations.
%
% In:
%   IndependentPath : platform-independent path; may contain forward and/or backward slashes
%                     (forward slashes generally preferred), and may refer to locations such as
%                     store:/ (the store path) or data:/ (one of the data paths); can also be a
%                     relative path
%
% Out:
%   SystemPath : system-specific path (with slashes corrected and locations resolved)
%                if multiple data paths are present, the one where the minimum number of
%                directories (and files) would have to be created to write to the given file is
%                selected.
%   
% Examples:
%   % translate a platform-independent reference to a subdirectory of the data path to one that is 
%   % recognized by the operating system (output might be, e.g., 'C:\Projects\mydata\test.mat')
%   env_translatepath('data:/projects/test.mat'); 
%
%   % resolve a reference to the storage directory
%   env_translatepath('store:/studyXY/result.mat')
%
%   % resolve a reference to the BCILAB root directory
%   env_translatepath('bcilab:/userscripts/myscript.m')
%
%   % resolve a reference to the user's home directory
%   env_translatepath('home:/myconfig.m')
%
%   % resolve a reference to the current temp directory (as specified in the startup options)
%   env_translatepath('temp:/output/001.mat')
%
%   % resolve a reference to the resources directory
%   env_translatepath('resources:/workspaces/testing.mat')
%
%
% See also:
%   env_startup, io_loadset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-29

global tracking;

% turn the path into a system-dependent one
filename = strrep(strrep(filename,'\',filesep),'/',filesep);

% resolve location references
if strncmp('store:',filename,6)
    filename = [tracking.paths.store_path filename(1+length('store:'):end)]; 
elseif strncmp('resources:',filename,10)
    filename = [tracking.paths.resource_path filename(1+length('resources:'):end)]; 
elseif strncmp('temp:',filename,5)
    filename = [tracking.paths.temp_path filename(1+length('temp:'):end)]; 
elseif strncmp('functions:',filename,10)
    filename = [tracking.paths.function_path filename(1+length('functions:'):end)]; 
elseif strncmp('bcilab:',filename,7)
    filename = [tracking.paths.bcilab_path filename(1+length('bcilab:'):end)]; 
elseif strncmp('dependencies:',filename,13)
    filename = [tracking.paths.dependency_path filename(1+length('dependencies:'):end)];
elseif strncmp('home:',filename,5)
    filename = [hlp_homedir filename(1+length('home:'):end)];
elseif strncmp('data:',filename,5)
    rest = filename(1+length('data:'):end);
    bestpath = 1; bestlen = -1;
    if length(tracking.paths.data_paths) > 1
        fpieces = hlp_split(rest,filesep);
        % find the data path that contains the longest prefix of the filename
        for pidx=1:length(tracking.paths.data_paths)
            p = tracking.paths.data_paths{pidx};
            % for each prefix of the filename (starting with the longest one)
            for k=length(fpieces):-1:0
                % check if the data path plus the first k pieces of the filename exists
                if exist([p sprintf([filesep '%s'],fpieces{1:k})],'file')
                    % found a match - check if it is a new length record among all our data paths...
                    if k>bestlen
                        bestlen = k;
                        bestpath = pidx;
                    end
                    break;
                end
            end
        end
    end
    % resolve the reference using that data path which matches most of the filename,
    % where, if multiple data paths are equally well suited, the first one of them is taken
    filename = [tracking.paths.data_paths{bestpath} rest];
end
