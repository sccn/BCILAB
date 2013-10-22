% set default config parameters

% ... additional parameters are permitted; see code/environment/env_startup.m for the possible
% options. If you run this bcilab distribution from multiple computers and need different settings
% for them, it is best to use conditional expressions like "if ispc", etc. if you share this bcilab
% distribution with other users, it is recommended to copy this file to your home directory (where
% it will be automatically found) or into a personal directory, and pass its file path as an
% argument to the "bcilab" function (e.g. cd /your/path/to/bcilab; bcilab
% /your/path/to/your/bcilab_config.m)

% this is your data path, i.e., where your studies and data sets are stored
data = {'/data/projects'};
data = [];

% this is a path where results can be stored (you must have write permission)
% (if empty, it is set to equal the data path)
store = [];

% this is the place where BCILAB will put its temporary files
% ideally, this is a fast local disk (e.g. an SSD disk), although network drives can be fast enough, 
% too (if this is set to [], bcilab does not cache its intermediate results on disk)
if exist('/data/cache','dir')    
    cache = {'dir','/data/cache','time',15};
elseif exist('/data/tmp','dir')
    cache = {'dir','/data/tmp','time',15};
elseif exist('/var/tmp','dir')
    cache = {'dir','/var/tmp','time',15};
else
    cache = [];
end
% comment out this line to enable caching of intermediate results
cache = [];

% This is the memory capacity that is being reserved for in-memory caching of intermediate results.
% For efficiency, it is important that this is at least large enough to hold a single data set.
% If this is smaller than one, it will be taken as the fraction of free system memory, otherwise it
% is in GB; A good value for workstations with 8 or more GB of memory is 0.25; for shared cluster 
% nodes with huge amounts of memory it is probably too generous (consider allocating a fixed amount).
% For under-powered laptops etc. a value of 0.5 allows for better memory utilization.
mem_capacity = 0.25;

% this is the cluster configuration; cluster computation is only partially tested so far, 
% especially on busy clusters where many things tend to go wrong.
% If you are starting your workers manually (in MATLAB: cd /path/to/bcilab; bcilab worker true),
% you can set this to something like {'engine','BLS','pool',{'hostname1:23547','hostname1:23548','hostname2:23547'}}
% if you are on a linux/UNIX machine talking to a linux cluster with MATLAB installed, 
% you can set it to something like {'engine','BLS','pool',{'hostname1','hostname2','hostname3'},'acquire',true},
% But to automatically acquire machines, you need to have a Linux .ssh identity set up so that 
% you can ssh into any of these machines without needing to type a password. Also, make sure that 
% these machines/ports are never accessible from the internet (except perhaps via SSH tunnels).
parallel = {'engine','local', 'pool', {'localhost:23547', 'localhost:23548', 'localhost:23549', 'localhost:23550', 'localhost:23551', 'localhost:23552', 'localhost:23553', 'localhost:23554'}};

% this controls how likely BCILAB is to cache things on disk; since disk access may be slower than
% recomputing the results again, this is the expected number of times you believe you'll need a
% given processed version of a data set again (e.g. a FIR highpassed version of a given data set).
data_reuses = 3;

% this is a directory for temporary results (e.g. IC decompositions); if this is empty, it will be
% set to a directory next to your BCILAB path
temp = [];

% whether to show the main BCILAB menu by default (if this is set to 'separate', the menu is 
% always detached from the EEGLAB main menu; otherwise it is under Tools>BCILAB)
menu = 'separate';

% whether to show experimental methods/options in the GUI
% (note that some of these may be unfinished features or otherwise very experimental)
show_experimental = false;

% whether to show guru-level options in the GUI by default (this can be toggled while editing)
show_guru = false;

% custom options to control how workers are acquired; parameters to par_getworkers_*
% these options are editable in the cluster settings dialog and are used when calling env_acquire_cluster (or clicking the "request cluster availability" button)
acquire_options = {'Hostnames',{'localhost'},'ShutdownTimeout',300};




