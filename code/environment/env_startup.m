function env_startup(varargin)
% Start the BCILAB toolbox, i.e. set up global data structures and load dependency toolboxes. 
% env_startup(Options...)
%
% Does all steps necessary for loading the toolbox -- the functions bcilab.m and eegplugin_bcilab.m
% are wrappers around this function which provide a higher level of convenience (configuration files
% in particular). Directly calling this function is not recommended.
%
% In:   
%  Options... : optional name-value pairs; allowed names are:
%
%               --- directory settings ---
%
%               'data':  Path where data sets are stored, used by data loading/saving routines.
%                        (default: path/to/bcilab/userdata)
%                        Note: this may also be a cell array of directories, in which case references 
%                              to data:/ are looked up in all of the specified directories, and the
%                              best match is taken.
%
%               'store': Path in which data shall be stored. Write permissions necessary (by default 
%                        identical to the data path)
%
%               'temp': temp directory (for misc outputs, e.g., AMICA models and dipole fits) 
%                       (default: path/to/bcilab-temp, or path/to/cache/bcilab_temp if a cache
%                       directory was specified)
%
%               'private': optional private plugin directory, separate from bcilab/*; can have
%                          the same directory structure as ~/.bcilab/ (default: [])
%
%               --- caching settings ---
%
%               'cache': Path where intermediate data sets are cached. Should be located on a fast 
%                        (local) drive with sufficient free capacity.
%                        * if this is left unspecified or empty, the cache is disabled
%                        * if this is a directory, it is used as the default cache location
%                        * a fine-grained cache setup can be defined by specifying a cell array of 
%                          cache locations, where each cache location is a cell array of name-value
%                          pairs, with possible names:
%                          'dir': directory of the cache location (e.g. '/tmp/bcilab_tmp/'), mandatory
%                          'time': only computations taking more than this many seconds may be stored 
%                                  in this location, but if a computation takes so long that another
%                                  cache location with a higher time applies, that other location is
%                                  preferred. For example, the /tmp directory may take computations
%                                  that take at least a minute, the home directory may take
%                                  computations that take at least an hour, the shared /data/results
%                                  location of the lab may take computations that take at least 12
%                                  hours (default: 30 seconds)
%                          'free': minimum amount of space to keep free on the given location, in GiB,
%                                   or, if smaller than 1, free is taken as the fraction of total
%                                   space to keep free. (default: 80)
%                          'tag': arbitrary identifier for the cache location (default: 'location_i', 
%                                 for the i'th location) must be a valid MATLAB struct field name,
%                                 only for display purposes
%
%               'mem_capacity': capacity of the memory cache (default: 2)
%                               if this is smaller than 1, it is taken as a fraction of the total
%                               free physical memory at startup time, otherwise it is in GB
%
%               'data_reuses' : estimated number of reuses of a data set being computed (default: 3)
%                               ... depending on disk access speeds, this determines whether it makes
%                               sense to cache the data set
%
%               --- parallel computing settings ---
%
%               'parallel' : parallelization options; cell array of name-value pairs, with names:
%                             'engine': parallelization engine to use, can be 'local', 
%                                       'ParallelComputingToolbox', or 'BLS' (BCILAB Scheduler) 
%                                       (default: 'local')
%                             'pool': node pool, cell array of 'host:port' strings; necessary for the 
%                                     BLS scheduler (default: {'localhost:23547','localhost:23548',
%                                     ..., 'localhost:23554'})
%                             'policy': scheduling policy function; necessary for the BLS scheduler 
%                                       (default: 'par_reschedule_policy')
%
%                             note: Parallel processing has so far received only relatively little
%                                   testing. Please use this facility at your own risk and report
%                                   any issues (e.g., starving jobs) that you may encounter.
%
%               'acquire_method' : name of the method to use to acquire worker processes (default:
%                                  'SSH')
%
%               'aquire_options' : Cell array of arguments as expected by par_getworkers_<Method>
%                                  (with bcilab-specific defaults for unspecified arguments)
%                                  (default: {})
%
%               'autokill_workers' : Whether to automatically kill workers when cluster resources
%                                    are released (default: true)
%
%               'worker' : whether this toolbox instance is started as a worker process or not; if 
%                          false, diary logging and GUI menus and popups are enabled. If given as a 
%                          cell array, the cell contents are passed on to the function par_worker,
%                          which lets the toolbox act as a commandline-less worker (waiting to 
%                          receive jobs over the network) (default: false)
%
%               --- misc settings ---
%
%               'menu' : create a menu bar (default: true) -- if this is set to 'separate', the BCILAB
%                        menu will be detached from the EEGLAB menu, even if run as a plugin.
%
%               'autocompile' : whether to try to auto-compile mex/java files (default: true)
%
%               'show_experimental' : whether to show methods and options marked 'experimental' in 
%                                     the GUIs (default: false)
%
%               'show_guru' : whether to show methods and options marked 'guru' in the GUIs
%                             (default: false)
%
%               'fingerprinting' : whether to calculate fingerprints of datasets; this allows to
%                                  manually curate dataset structs "EEGLAB-style" prior to using
%                                  them in offline analysis functions (e.g., bci_train). If
%                                  disabled, all dataset processing must happen in custom filters,
%                                  loaders, or must have been committed into .set files. Disabling
%                                  this can give a speed advantage on extremely large data sets.
%                                  (default: true)
%
%               'disabled_caches' : cell array of disabled disk caches; see hlp_diskcache lines at
%                                   the bottom of this file for available caches (default: {})
%
%
% Examples:
%  Note that env_startup is usually not being called directly; instead the bcilab.m function in the
%  bcilab root directory forwards its arguments (and variables declared in a config script) to this 
%  function.
%
%  % start BCILAB with a custom data directory, and a custom cache directory
%  env_startup('data','C:\Data', 'cache','C:\Data\Temp');
%
%  % as before, but specify multiple data paths that are being fused into a common directory view
%  % where possible (in case of ambiguities, the earlier directories take precedence)
%  env_startup('data',{'C:\Data','F:\Data2'}, 'cache','C:\Data\Temp');
%
%  % start BCILAB with a custom data and storage directory, and specify a cache location with some
%  % additional meta-data (namely: only cache there if a computation takes at least 60 seconds, and 
%  % reserve 15GB free space
%  env_startup('data','/media/data', 'store','/media/data/results', 'cache',{{'dir','/tmp/tracking','time',60,'free',15}});
% 
%  % as before, but make sure that the free space does not fall below 20% of the disk
%  env_startup('data','/media/data', 'store','/media/data/results', 'cache',{{'dir','/tmp/tracking','time',60,'free',0.2}});
%
%  % start BCILAB and set up a very big in-memory cache for working with large data sets (of 16GB)
%  env_startup('mem_capacity',16)
%
%  % start BCILAB but prevent the main menu from popping up
%  env_startup('menu',false)
%
%  % start BCILAB and specify which parallel computation resources to use; this assumes that the 
%  % respective hostnames are reachable from this computer, and are running MATLAB sessions which 
%  % execute a command similar to: cd /your/path/to/bcilab; bcilab('worker',true); par_worker;
%  env_startup('parallel',{'engine','BLS', 'pool',{'computer1','computer2','computer3'}}
%  
%
%  % start the toolbox as a worker
%  env_startup('worker',true)
%
%  % start the toolbox as worker, and pass some arguments to par_worker (making it listen on port
%  % 15456, using a portrange of 0, and using some custom update-checking arguments
%  env_startup('worker',{15456,0,'update_check',{'/data/bcilab-0.9-beta2b/build/bcilab','/mnt/remote/bcilab-build/bcilab'}})
%
% See also:
%  env_load_dependencies, env_translatepath
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% benchmarking
t_prelaunch = tic;

% determine the BCILAB core directories
if ~isdeployed
    tmpdir = path_normalize(fileparts(mfilename('fullpath')));
    delims = strfind(tmpdir,filesep);
    base_dir = tmpdir(1:delims(end-1));
else
    % in deployed mode, we walk up until we find the BCILAB base directory
    tmpdir = pwd;
    delims = strfind(tmpdir,filesep);
    disp(['Launching from directory ' tmpdir ' ...']);
    for k=length(delims):-1:1
        base_dir = tmpdir(1:delims(k));        
        if exist([base_dir filesep 'code'],'dir')
            success = true; %#ok<NASGU>
            break; 
        end
    end
    if ~exist('success','var')
        error('Could not find the ''code'' directory; make sure that the binary is in a sub-directory of a full BCILAB distribution.'); end
end

function_dir = [base_dir 'code'];
dependency_dir = [base_dir 'dependencies'];
resource_dir = [base_dir 'resources'];
script_dir = [base_dir 'userscripts'];
build_dir = [base_dir 'build'];

        
% add them all to the MATLAB path (except for dependencies, which are loaded separately)
if ~isdeployed
    
    % remove existing BCILAB path references
    ea = which('env_add');
    if ~isempty(ea)
        % get the bcilab root directory that's currently in the path
        bad_path = ea(1:strfind(ea,'dependencies')-2);
        % remove all references
        paths = strsplit(path,pathsep);
        retain = cellfun('isempty',strfind(paths,bad_path));
        path(sprintf(['%s' pathsep],paths{retain}));
        if ~all(retain)
            disp('BCILAB sub-directories have been detected in the MATLAB path, removing them.'); end
    end

    % add core function paths
    addpath(genpath(function_dir));
    if exist(build_dir,'dir')
        addpath(build_dir); end
    evalc('addpath(genpath(script_dir))');
    evalc('addpath([base_dir ''userdata''])');
    
    % remove existing eeglab path references, if BCILAB is not itself contained as a plugin in this 
    % EEGLAB distribution
    ep = which('eeglab');
    if ~isempty(ep) && isempty(strfind(mfilename('fullpath'),fileparts(which('eeglab'))))
        paths = strsplit(path,pathsep);
        ep = strsplit(ep,filesep);
        retain = cellfun('isempty',strfind(paths,ep{end-1}));
        path(sprintf(['%s' pathsep],paths{retain}));
        if ~all(retain)
            disp('  The previously loaded EEGLAB path has been replaced.'); end
    end
end

if hlp_matlab_version < 706
    disp('Note: Your version of MATLAB is not supported by BCILAB any more. You may try BCILAB version 0.9, which supports old MATLAB''s back to version 2006a.'); end    

% get options
opts = hlp_varargin2struct(varargin,'data',[],'store',[],'cache',[],'temp',[],'private',[],'mem_capacity',2,'data_reuses',3, ...
    'parallel',{'use','local'}, 'menu',true, 'configscript','', 'worker',false, 'autocompile',true, 'acquire_method','SSH', 'acquire_options',{}, ...
    'autokill_workers',true,'show_experimental',false,'show_guru',false,'fingerprinting',true,'disabled_caches',{});

% load all dependencies, recursively...
disp('Loading BCILAB dependencies...');
env_load_dependencies(dependency_dir,opts.autocompile);
if exist(env_translatepath('home:/.bcilab/code/dependencies'),'dir')
    env_load_dependencies(env_translatepath('home:/.bcilab/code/dependencies'),opts.autocompile); end
if ~isempty(opts.private) && exist(env_translatepath([opts.private '/code/dependencies']),'dir')
    env_load_dependencies(env_translatepath([opts.private '/code/dependencies']),opts.autocompile); end

if ischar(opts.worker)
    try
        disp(['Evaluating worker argument: ' opts.worker]);
        opts.worker = eval(opts.worker); 
    catch
        disp('Failed to evaluate worker; setting it to false.');
        opts.worker = false; 
    end
elseif ~isequal(opts.worker,false)
    fprintf('Worker was given as a %s with value %s\n',class(opts.worker),hlp_tostring(opts.worker));
end
if ~isequal(opts.worker,false) && ~iscell(opts.worker)
    error('The given worker argument must be a cell array or false.'); end

if ischar(opts.parallel)
    try
        disp(['Evaluating parallel argument: ' opts.parallel]);
        opts.parallel = eval(opts.parallel); 
    catch
        disp('Failed to evaluate worker; setting it to empty.');
        opts.parallel = {};
    end
end

% process data directories
if isempty(opts.data)
    opts.data = {}; end
if ~iscell(opts.data)
    opts.data = {opts.data}; end
for d = 1:length(opts.data)
    opts.data{d} = path_normalize(opts.data{d}); end
if isempty(opts.data)
    opts.data = {[base_dir 'userdata']}; end

% process store directory
if isempty(opts.store)
    opts.store = opts.data{1}; end
opts.store = path_normalize(opts.store);

% process cache directories
if isempty(opts.cache)
    opts.cache = {}; end
if ischar(opts.cache)
    opts.cache = {{'dir',opts.cache}}; end
if iscell(opts.cache) && ~isempty(opts.cache) && ~iscell(opts.cache{1})
    opts.cache = {opts.cache}; end
for d=1:length(opts.cache)
    opts.cache{d} = hlp_varargin2struct(opts.cache{d},'dir','','tag',['location_' num2str(d)],'time',30,'free',80); end
% remove entries with empty dir
opts.cache = opts.cache(cellfun(@(e)~isempty(e.dir),opts.cache));
for d=1:length(opts.cache)
    % make sure that the BCILAB cache is in its own proper sub-directory
    opts.cache{d}.dir = [path_normalize(opts.cache{d}.dir) filesep 'bcilab_cache'];
    % create the directory if necessary
    if ~isempty(opts.cache{d}.dir) && ~exist(opts.cache{d}.dir,'dir')
        try
            io_mkdirs([opts.cache{d}.dir filesep],{'+w','a'});
        catch
            disp(['cache directory ' opts.cache{d}.dir ' does not exist and could not be created']);
        end
    end
end

% process temp directory
if isempty(opts.temp)    
    if ~isempty(opts.cache)
        opts.temp = [fileparts(opts.cache{1}.dir) filesep 'bcilab_temp']; 
    else
        opts.temp = [base_dir(1:end-1) '-temp'];
    end
end
opts.temp = path_normalize(opts.temp);
try
    io_mkdirs([opts.temp filesep],{'+w','a'});
catch
    disp(['temp directory ' opts.temp ' does not exist and could not be created.']);
end

% set global variables
global tracking
tracking.paths = struct('bcilab_path',{base_dir(1:end-1)}, 'function_path',{function_dir}, 'data_paths',{opts.data}, 'store_path',{opts.store}, 'dependency_path',{dependency_dir},'resource_path',{resource_dir},'temp_path',{opts.temp}, 'private_path',{opts.private});
for d=1:length(opts.cache)
    location = rmfield(opts.cache{d},'tag');
    % convert GiB to bytes
    if location.free >= 1 
        location.free = location.free*1024*1024*1024; end    
    try
        warning off MATLAB:DELETE:Permission;
        % probe the cache locations...
        import java.io.*; 
        % try to add a free space checker (Java File object), which we use to check the quota, etc.
        location.space_checker = File(opts.cache{d}.dir);
        [pid,pname] = hlp_processid;
        filename = [opts.cache{d}.dir filesep '__probe_cache_ ' num2str(pid) '_' pname '__.sto'];
        if exist(filename,'file') 
            delete(filename); end
        oldvalue = location.space_checker.getFreeSpace;
        testdata = double(rand(1024)); %#ok<NASGU>
        objinfo = whos('testdata');
        % do a quick read/write test
        t0=tic; io_save(filename,'testdata'); location.writestats = struct('size',{0 objinfo.bytes},'time',{0 toc(t0)});
        t0=tic; io_load(filename); location.readstats = struct('size',{0 objinfo.bytes},'time',{0 toc(t0)});
        newvalue = location.space_checker.getFreeSpace;
        if exist(filename,'file') 
            delete(filename); end
        % test if the space checker works, and also get some quick measurements of disk read/write speeds
        if newvalue >= oldvalue
            location = rmfield(location,'space_checker'); end        
        % and turn the free space ratio into an absolute value
        if location.free < 1
            location.free = location.free*location.space_checker.getTotalSpace; end
    catch e
        disp(['Could not probe cache file system speed; reason: ' e.message]);
        location.writestats = struct('size',{0 1024},'time',{0 0.01});
        location.readstats = struct('size',{0 1024},'time',{0 0.01});
        if location.free < 1
            location.free = 60*2^20; end % if space check fails, assume that 60GB are free
    end
    tracking.cache.disk_paths.(opts.cache{d}.tag) = location;     
end
if opts.mem_capacity < 1
    free_mem = hlp_memfree();
    tracking.cache.capacity = round(opts.mem_capacity * free_mem);
    if free_mem < 1024*1024*1024
        sprintf('Warning: You have less than 1 GB of free memory (reserving %.0f%% = %.0fMB for data caches).\n',100*opts.mem_capacity,tracking.cache.capacity/(1024*1024));
        sprintf('         This will severely impact the offline processing speed of BCILAB.\n');
        sprintf('         You may force a fixed amount of cache cacpacity by assinging a value greater than 1 (in GB) to the ''mem_capacity'' variable in your bcilab_config.m.'); 
    end
else
    tracking.cache.capacity = opts.mem_capacity*1024*1024*1024;
end
tracking.cache.reuses = opts.data_reuses;
tracking.cache.data = struct();
tracking.cache.sizes = struct();
tracking.cache.times = struct();
if ~isfield(tracking.cache,'disk_paths')
    tracking.cache.disk_paths = struct(); end
% initialize stack mechanisms
tracking.stack.base = struct('disable_expressions',false);
% set parallelization settings
tracking.parallel = hlp_varargin2struct(opts.parallel, ...
    'engine','local', ...
    'pool',{'localhost:23547','localhost:23548','localhost:23549','localhost:23550','localhost:23551','localhost:23552','localhost:23553','localhost:23554'}, ...
    'policy','par_reschedule_policy',...
    'verbosity',0);
tracking.acquire_method = opts.acquire_method;
tracking.acquire_options = opts.acquire_options;
tracking.autokill_workers = opts.autokill_workers;
tracking.configscript = opts.configscript;
% set GUI settings
tracking.gui.show_experimental = opts.show_experimental;
tracking.gui.show_guru = opts.show_guru;
try
    cd(script_dir);
catch
end

% set up some microcache domain properties
hlp_microcache('arg','lambda_equality','proper');
hlp_microcache('spec','group_size',5);
hlp_microcache('findfunction','lambda_equality','fast','group_size',5);
hlp_microcache('verybig','max_key_size',2^30,'max_result_size',2^30);

% set up some fine-grained disk-cache locations
% small designed filter models in signal processing functions (crucial, especially due to toolbox licenses; frequently reused)
hlp_diskcache('filterdesign','folder',opts.temp,'subdir','filterdesign','exactmatch_cutoff',0,'bypass',ismember('filterdesign',opts.disabled_caches));
% ICA solutions; small and expensive to recompute
hlp_diskcache('icaweights','folder',opts.temp,'subdir','icaweights','exactmatch_cutoff',0,'spot_hashing',true,'bypass',ismember('icaweights',opts.disabled_caches));
% dipole fits: small and expensive to recompute
hlp_diskcache('dipfits','folder',opts.temp,'subdir','dipfits','exactmatch_cutoff',0,'bypass',ismember('dipfits',opts.disabled_caches));
% feature-extraction models: fairly small, but produced in large quantities and fast to recompute
hlp_diskcache('featuremodels','folder',opts.temp,'subdir','featuremodels','exactmatch_cutoff',0,'spot_hashing',true,'bypass',ismember('featuremodels',opts.disabled_caches));
% machine learning models: can be large, can be expensive to recompute, but produced in large quantities
hlp_diskcache('predictivemodels','folder',opts.temp,'subdir','predictivemodels','exactmatch_cutoff',0,'spot_hashing',true,'bypass',ismember('predictivemodels',opts.disabled_caches));
% miscellaneous statistics; currently rarely used
hlp_diskcache('statistics','folder',opts.temp,'subdir','statistics','exactmatch_cutoff',0,'bypass',ismember('statistics',opts.disabled_caches));
% montage-related data: small but expensive to recompute (frequently reused)
hlp_diskcache('montages','folder',opts.temp,'subdir','montages','exactmatch_cutoff',0,'bypass',ismember('montages',opts.disabled_caches));
% more montage-related data: small but expensive to recompute (frequently reused)
hlp_diskcache('montage_quality','folder',opts.temp,'subdir','montages','exactmatch_cutoff',0,'spot_hashing',true,'bypass',ismember('montage_quality',opts.disabled_caches));
% intermediate results: produced in huge quantities
hlp_diskcache('intermediate','folder',opts.temp,'subdir','intermediate','bypass',ismember('intermediate',opts.disabled_caches));
% feature representations: used by some expensive methods; moderately expensive to recompute per byte, produced in huge quantities
hlp_diskcache('features','folder',opts.temp,'subdir','features','exactmatch_cutoff',0,'spot_hashing',true,'bypass',ismember('features',opts.disabled_caches));
% generic data: currently rarely used
hlp_diskcache('general','folder',opts.temp,'subdir','general','bypass',ismember('general',opts.disabled_caches));
% fine-grained data: currently rarely used
hlp_diskcache('finegrained','folder',opts.temp,'subdir','finegrained','bypass',ismember('finegrained',opts.disabled_caches));

% set up fingerprinting options
if ~opts.fingerprinting
    exp_eval(exp_set('fingerprint_create',false));
    exp_eval(exp_set('fingerprint_check',false));
end


% show toolbox status
fprintf('\n');
if ~isempty(opts.disabled_caches)
    fprintf('The following caches are disabled: %s.\n',hlp_tostring(opts.disabled_caches)); end
if ~opts.fingerprinting
    fprintf('Caution: data fingerprinting has been disabled. You cannot pass manually edited datasets into offline analysis functions.\n'); end
fprintf('\n');
disp(['code is in ' function_dir]);
datalocs = [];
for d = opts.data
    datalocs = [datalocs d{1} ', ']; end %#ok<AGROW>
disp(['data is in ' datalocs(1:end-2)]);
disp(['results are in ' opts.store]);
if ~isempty(opts.cache)
    fnames = fieldnames(tracking.cache.disk_paths);
    for f = 1:length(fnames)
        if f == 1
            disp(['cache is in ' tracking.cache.disk_paths.(fnames{f}).dir ' (' fnames{f} ')']);
        else
            disp(['            ' tracking.cache.disk_paths.(fnames{f}).dir ' (' fnames{f} ')']);
        end
    end    
else
    disp('cache is disabled');
end
disp(['temp is in ' opts.temp]);
if ~isempty(opts.private)
    disp(['private plugins are in ' opts.private]); end
fprintf('\n');


% turn off a few warnings
warning off MATLAB:structOnObject
warning off MATLAB:log:logOfZero
warning off MATLAB:divideByZero %#ok<RMWRN>
warning off MATLAB:RandStream:ReadingInactiveLegacyGeneratorState % for GMMs....

if isequal(opts.worker,false) || isequal(opts.worker,0)
    % --- regular mode ---   
    
    try 
        % create directories in the user's .bcilab folder...
        home_basedir = [hlp_homedir filesep '.bcilab' filesep];
        home_codedirs = {['code' filesep 'filters'],['code' filesep 'dataset_editing'], ...
            ['code' filesep 'machine_learning'], ['code' filesep 'paradigms'], ['code' filesep 'scripts']};
        home_miscdirs = {'models','approaches','code',['code' filesep 'dependencies'],'logs',['logs' filesep 'workers']};
        for d = [home_codedirs home_miscdirs]
            try
                subdir = [home_basedir d{1}];
                if ~exist(subdir,'dir')
                    mkdir(subdir); end
            catch e
                disp(['Could not create directory: ' subdir ': ' e.message]);
            end
        end
        % and add the code directories to the path
        for d = home_codedirs
            addpath(genpath([home_basedir d{1}])); end
        % add the env_add.m file to the dependencies folder
        add_file = [home_basedir 'code' filesep 'dependencies' filesep 'env_add.m'];
        if ~exist(add_file,'file')
            try                
                f = fopen(add_file,'w');
                fwrite(f,' ');
                fclose(f);
            catch e
                disp(['Could not create file ' add_file ': ' e.message]);
            end
        end
    catch e
        disp(['Error trying to set up directories in ~/.bcilab: ' e.message]);
    end

    % set up logfile
    if ~exist([hlp_homedir filesep '.bcilab'],'dir')
        if ~mkdir(hlp_homedir,'.bcilab');
            disp('Cannot create directory .bcilab in your home folder.'); end
    end
    tracking.logfile = env_translatepath(sprintf('home:/.bcilab/logs/bcilab_%s-%s.log',hlp_hostname,strrep(datestr(now),':','.')));
    try 
        if exist(tracking.logfile,'file')
            warning off MATLAB:DELETE:Permission
            delete(tracking.logfile); 
            warning on MATLAB:DELETE:Permission
        end
        diary(tracking.logfile); 
    catch e
        disp(['Error trying to set up the logfile: ' e.message]);
    end   
    
    % add private plugin directories to the path
    if ~isempty(opts.private)
        try
            private_codedirs = {['code' filesep 'filters'],['code' filesep 'dataset_editing'], ...
                ['code' filesep 'machine_learning'], ['code' filesep 'paradigms'], ['code' filesep 'scripts']};
            % and add the code directories to the path
            for d = private_codedirs
                tmpdir  = env_translatepath([opts.private filesep d{1}]);
                if exist(tmpdir,'dir')
                    addpath(genpath(tmpdir)); end
            end
        catch e
            disp(['Could not add private code directories to path: ' e.message]);
        end
    end
    
    % create a menu
    if ~(isequal(opts.menu,false) || isequal(opts.menu,0))
        try
            env_showmenu('forcenew',strcmp(opts.menu,'separate'));
        catch e
            disp('Could not open the BCILAB menu; traceback: ');
            env_handleerror(e);
        end
    end
    
    % display a version reminder
    bpath = hlp_split(env_translatepath('bcilab:/'),filesep);
    if ~isempty(strfind(bpath{end},'-stable'))
        if isdeployed
            disp(' This is the stable version.');
        else
            disp(' This is the stable version - please keep this in mind when editing.');
        end
    elseif ~isempty(strfind(bpath{end},'-devel'))
        if isdeployed
            disp(' This is the DEVELOPER version.');
        else
            try
                cprintf([0 0.4 1],'This is the DEVELOPER version.\n');
            catch
                disp(' This is the DEVELOPER version.');
            end
        end
    end

    disp([' Welcome to the BCILAB toolbox on ' hlp_hostname '!'])
    fprintf(' Launch time was: %.1f seconds\n',toc(t_prelaunch));
    fprintf('\n');

else
    disp('Now entering worker mode...');
    fprintf(' Launch time was: %.1f seconds\n',toc(t_prelaunch));
    fprintf('\n');
    
    % -- worker mode ---
    % disable standard dialogs in the workers
    if ~isdeployed
        addpath(env_translatepath('dependencies:/disabled_dialogs')); end
    
    % close EEGLAB main menu
    mainmenu = findobj('Tag','EEGLAB');
    if ~isempty(mainmenu)
        close(mainmenu); end
    drawnow;
    
    % translate the options
    if isequal(opts.worker,true)
        opts.worker = {}; end
    if ~iscell(opts.worker)
        opts.worker = {opts.worker}; end
    % start!
    par_worker(opts.worker{:});
end

try
    % pretend to invoke the dependency list so that the compiler finds it...
    dependency_list;
catch
end
    
    
% normalize a directory path
function dir = path_normalize(dir)
dir = strrep(strrep(dir,'\',filesep),'/',filesep);
if dir(end) == filesep
   dir = dir(1:end-1); end


% Split a string according to some delimiter(s). % Not as fast as hlp_split (and doesn't fuse 
% delimiters), but works without bsxfun.
function strings = strsplit(string, splitter)
ix = strfind(string, splitter);
strings = cell(1,numel(ix)+1);
ix = [0 ix numel(string)+1];
for k = 2 : numel(ix)
    strings{k-1} = string(ix(k-1)+1:ix(k)-1); end
