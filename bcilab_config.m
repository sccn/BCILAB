% BCILAB configuration file
% sets default config parameters

% Additional parameters than those listed here are permitted; see code/environment/env_startup.m for
% the possible options. If you run this bcilab distribution from multiple computers and need
% different settings for them, it is possible to use conditional expressions like "if ispc", etc. if
% you share this bcilab distribution with other users, it is recommended to copy this file to your
% home directory (under .bcilab, where it will be automatically found) or into a personal directory,
% and pass its file path as an argument to the "bcilab" function (e.g. cd /your/path/to/bcilab;
% bcilab /your/path/to/your/bcilab_config.m)

% This file is loaded by bcilab.m -- note that you can have a version of the file in your ~/.bcilab/
% directory, which will take precedence over a file of same name placed in the bcilab folder (see
% bcilab.m), and that you can specify a custom config file to load. You can determine your home
% directory by calling hlp_homedir after bcilab has loaded.
%
% The options in this configuration file can also be edited via the GUI (under Settings); the GUI
% always edits the file that was used when the toolbox was loaded.

% This is your data path, i.e., where your studies and data sets are stored; you can also have multiple
% possible paths in a cell array (if empty, will default to your bcilab:/userdata folder).
data = []; % e.g., {'/data/myprojects','bcilab:/userdata','home:/mydata'};

% This is a path where results can be stored (you must have write permission).
% (if empty, it is set to equal the data path)
store = [];

% This is the place where BCILAB will cache intermediate results (mostly data sets). Ideally, this
% is a fast local disk (e.g. an SSD disk), although network drives can be fast enough, too (if this
% is set to [], bcilab does not cache its intermediate results on disk).
if exist('/data/cache','dir')    
    cache = {'dir','/data/cache','time',15};
elseif exist('/data/tmp','dir')
    cache = {'dir','/data/tmp','time',15};
elseif exist('/var/tmp','dir')
    cache = {'dir','/var/tmp','time',15};
else
    cache = [];
end
% Comment out this line to enable caching of intermediate results in one of the above default linux
% temp folders.
cache = [];

% This is the memory capacity that is being reserved for in-memory caching of intermediate results.
% For efficiency, it is important that this is at least large enough to hold a single data set.
% If this is smaller than one, it will be taken as the fraction of free system memory, otherwise it
% is in GB; A good value for workstations with 8 or more GB of memory is 0.25; for shared cluster 
% nodes with huge amounts of memory it is probably too generous (consider allocating a fixed amount).
% For under-powered laptops etc. a value of 0.5 allows for better memory utilization.
mem_capacity = 0.25;

% This is the cluster configuration; cluster computation is only partially tested so far, especially
% on busy clusters where machines go down or become unavailable for various reasons. If you are
% starting your workers manually (in MATLAB: cd /path/to/bcilab; bcilab worker true), you can set
% this to something like {'engine','BLS','pool',{'myhost1:23547','myhost1:23548','myhost2:23547'}}
% where the hostnames/ports are those reported by the respective workers. If you are on a linux/UNIX
% machine talking to a linux cluster with MATLAB installed, it is usually more convenient to instead
% set the acquire_options (best done in the GUI's cluster settings) and acquire workers by the click
% of the acquire-cluster jbutton in the GUI (which sets the parallel setting based on which machines
% are available). But to automatically acquire machines, you need to have a Linux .ssh identity set
% up so that you can ssh into any of these machines without needing to type a password. Note that
% machines/ports that run workers should never be openly accessible from the internet or unsafe
% networks since it is possible to remote-control the machine this way.
parallel = {'engine','local', 'pool', {'localhost:23547', 'localhost:23548', 'localhost:23549', 'localhost:23550', 'localhost:23551', 'localhost:23552', 'localhost:23553', 'localhost:23554'}};

% This determines how likely BCILAB is to cache data sets on disk; since disk access may be slower
% than recomputing the results again, this is the expected number of times you believe you'll need a
% given processed version of a data set again (e.g. a high-pass filtered version of a given data
% set).
data_reuses = 3;

% This is a directory for temporary results; if this is empty, it will be set to a directory that is
% your BCILAB installation path with -temp appended.
temp = [];

% If you have private plugins that you manage separately from the bcilab directory tree, you can set
% the path to these here (using the same sub-directory structure as ~/.bcilab/, see that folder for
% reference).
private = [];

% Whether to show the main BCILAB menu by default (if this is set to 'separate', the menu is 
% always detached from the EEGLAB main menu; otherwise it is under Tools>BCILAB).
menu = 'separate';

% Whether to show experimental methods/options in the GUI(note that some of these may be unfinished
% features or otherwise very experimental).
show_experimental = true;

% Whether to show guru-level options in the GUI by default (this can be toggled while editing).
% Unless you are an experienced engineer it is recommended to keep these options hidden to unclutter
% the options menus.
show_guru = true;

% Custom options to control how workers are acquired; parameters to par_getworkers_* these options
% are editable in the cluster settings dialog and are used when calling env_acquire_cluster (or
% clicking the "request cluster availability" button).
acquire_method = 'SSH';
acquire_options = {'Hostnames',{'localhost'},'ShutdownTimeout',300};

% Whether to automatically kill workers when cluster resources are being released
% Note that workers typically shut themselves down within a preset idle time after having been released,
% so this option only serves to expedite the process (e.g., to start with fresh workers after a code change)
autokill_workers = true;


