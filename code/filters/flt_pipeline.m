function signal = flt_pipeline(varargin)
% Configurable preprocessing pipeline for most BCI paradigms.
% Signal = flt_pipeline(Signal, Stages...)
%
% Most BCI paradigms contain a sequence of signal processing steps, many of which are common to a
% variety of paradigms (for example, spatial and spectral filtering). The default Signal Processing
% pipeline allows to restrict the raw sensor signals to the contents of interest, such as, for
% example, specific frequency bands (e.g., the alpha band), or specific spatial areas (e.g. signals
% emitted from the motor cortex). This can be viewed as a way to improve the signal-to-noise ratio
% (for some notion of what is the signal) of the data, or as a way to express and incorporate prior
% knowledge into the paradigm.
%
%
% Direct Use
% ==========
%
% With no options, flt_pipeline returns the input signal unprocessed. By giving options, various
% pipeline stages can be enabled. By default, the stages are automatically ordered according to 
% whatever constraints and ordering hints they declare (which covers 99% of cases, though the order
% can be manually and selectively overridden), and the order in which they are passed to flt_pipeline
% does not matter. If a pipeline stage is listed multiple times, only the last assignment is 
% effective.
%
% Each pipeline stage can be referred to by various names (among other for backwards compatibility),
% including a) the name of the function (without the flt_ / set_ prefix), b) the name that shows 
% up in the GUI (e.g. SpectralTransform for flt_fourier), and c) in some cases a legacy short-cut 
% name, e.g. 'epoch' for set_makepos).
%
% The canonical way of enabling pipeline stages is by passing the name of the stage, followed by a
% cell array containing the list of arguments to be passed to the respective function, excluding the
% data set itself (which would come from whatever stage was before). Thus, the argument format is
% determined by the respective stage function (see examples in the function's help). An example is
% flt_pipeline('signal',eeg, 'iir',{'Frequencies',[0.5 2],'Mode','highpass'})
% 
% In addition, there are some shortcuts if only a single parameter that is not a cell array is
% passed to the respective stage. In this case, the parameter can usually be passed without the
% enclosing cell array, for example as in flt_pipeline('signal',eeg,'resample',300,'iir',[5 7 25 30]). 
% It is also possible to turn on or off a pipeline stage by passing the string 'on' or 'off'.
%
% The function automatically scans the filter and dataset_editing folders and passes the arguments
% on to the respective functions.
%
%
% Use in a BCI paradigm function
% ==============================
%
% As many BCI paradigms make use of a signal processing pipeline, flt_pipeline is most frequently
% used within BCI paradigms (within their preprocessing code). BCI paradigms might expose a subset
% of the filter pipeline's arguments as their own user parameters, or they might expose the entire
% parameter set of flt_pipeline to the user, typically with some pre-defined paradigm-specific 
% settings, so that the user can customize the entire pipeline. In particular, the function 
% para_dataflow is a template BCI paradigm which exposes the entire pipeline as its 'flt' or 
% 'SignalProcessing' argument (this is a cell array of flt_pipeline arguments), and the majority of 
% BCI paradigms are wrappers around this function (only changing some default stages).
%
% In this case, there is already a default assignment for a subset of the pipeline stage parameters,
% and the user can selectively override the defaults, for example, by passing 
% 'SignalProcessing',{'EpochExtraction',{...}} to change only the epoch extraction behavior of the
% paradigm. To disable a stage that is by default enabled by a paradigm, simply pass [],
% e.g. 'SignalProcessing',{'Resampling',[]} or 'off' to disable the respective stage. To turn a 
% stage on without passing any particular parameters, pass either {} or 'on'.
%
%
% Writing compatible stages
% =========================
%
% To be a compatible stage recognized by flt_pipeline, a function needs to:
% a) be in the directory code/filters or code/dataset_editing and begin with flt_ or set_
% b) use the argument declaration system (i.e., arg_define) to declare its arguments
% c) declare at least one argument named 'signal', which is the signal to be processed
%    and when returning a second output (the state), accept an optional argument called 'state'
% d) preferably use the expression system (i.e. exp_beginfun and exp_endfun)
% e) preferably declare at least a minimal set of ordering hints (using declare_properties)
% f) preferably declare a human-readable 'name' property in CamelCase (using declare_properties)
%
% 
% Filter Ordering
% ===============
%
% Each filter in BCILAB declares a small set of ordering hints which
% describes its required or preferred ordering with respect to other
% filters. This is to save users the hassle of figuring out the optimal
% filter order themselves in the majority of cases.
%
% For example, a filter that can only be applied on continuous but
% not segmented data would specify that it cannot follow the filter stage
% that turns a continuous signal into a segmented signal (which is the 
% function set_makepos). Generally, these ordering hints are expressed
% in the filter function (e.g. flt_iir) using a line as in the example:
%
%   declare_properties('cannot_follow','set_makepos', 'follows','flt_resample')
%  
% Aside from the two hard ordering constraints 'cannot_follow' and
% 'cannot_precede', a third one is the hard 'depends' constraint, which 
% expresses that a particular filter stage must have been applied
% beforehand -- for example, dipole fitting (for independent components) 
% can only be performed after independent components have been derived in 
% the first place, which is done by the filter flt_ica.
%
% Lastly, a filter may express also "soft" ordering preferences (and many
% do), for example because certain orders may be computationally more 
% efficient or because a higher-quality result can be attained in a certain 
% order. The two soft preferences are 'precedes' and 'follows'. For
% example, many filters prefer to be applied after the resampling stage
% (flt_resample) rather than before, in order to save computation time 
% (as resampling usually reduces the amount of data to process). Another 
% example are frequency-domain transforms (flt_fourier and flt_coherence), 
% which are most useful when applied to data that is spatially filtered
% rather than to the raw channel signals. Thus, these filters prefer to be
% applied after at least flt_ica and flt_project. Because many filters
% contribute such constraints or preferences (which can be visualized as 
% directed edges between filter nodes), the final order is fairly well 
% determined. 
%
% If a different ordering is required than the default (which can be looked
% up in the GUI edit panel of any paradigm that uses flt_pipeline -- that is
% most of them), then one many partially override the order by passing a 
% cell array of filter names as the parameter FilterOrdering. This is
% interpreted as a set of additional hard ordering constraints between any
% pair of filters in the list. Generally, the user-specified ordering and
% the hard constraints take precedence over the preferences, but hard
% constraints can not conflict (e.g., one cannot force a continuous-data
% filter to be applied after set_makepos but gets an error instead). 
%
% Lastly, only filters that specify at least one ordering hint are managed
% by flt_pipeline (possibly a reason why a hastily written filter may not
% show up in the GUI).
%
%
% In:
%   Signal     : a data set to be processed (as, e.g., loaded with io_loadset)
%
%   Stages... : optional name-value pairs that enable various stages of the default preprocessing
%               pipeline. In the following, the most common default stages are listed with some 
%               examples:
%
%                'Resampling'/'resample'/'srate' (see flt_resample):
%                   resample to the given sampling rate, in Hz; 
%
%                   Example: resample to 200 Hz
%                   flt_pipeline(X,...,'resample',200,...)
%
%
%                'Rereferencing'/'reref'/'ref' (see flt_reref):
%                   re-reference the data to a set of channel(s); 
%
%                   Example: resample to the average of 'TP7' and 'TP8' channels
%                   flt_pipeline(X,...,'Rereferencing',{{'TP7','TP8'}},...)
%
%                   Example: do a common average reference
%                   flt_pipeline(X,...,'Rereferencing',{[]},...)
%
%
%                'ICA'/'ica' (see flt_ica):
%                   annonatate the data set with an independent component analysis decomposition
%
%                   Example: use the default settings
%                   flt_pipeline(X,...,'ICA','on',...), 
%
%                   Example: use a 3-models amica, and use 16 nodes on the cluster, using 'hardcore' cleaning
%                   flt_pipeline(X,...,'ICA',{'Variant',{'amica','NumModels',3,'NumProcessors',16},'CleaningLevel','hardcore'}...)
%
%
%                'ChannelSelection'/'selchans'/'channels' (see flt_selchans): 
%                   select a channel subset, typically a cell-string array 
%                   
%                   Example: select channels C3 and C4
%                   flt_pipeline(X,...,'ChannelSelection',{{'C3' 'C4'}},...)
%
%
%                'SurfaceLaplacian'/'laplace' (see flt_laplace):
%                   simple Hjorth-style surface laplacian
%
%                   Example: use the default settings
%                   flt_pipeline(X,...,'laplace','on',...)
%
%
%                'IIRFilter'/'iir' (see flt_iir): 
%                   apply an IIR-based frequency filter
%
%                   Example: implement an 8-30 Hz band-pass filter
%                   flt_pipeline(X,...,'iir',[7 8 29 31],...)
%
%                   Example: implement an 1 Hz high-pass filter (with generous transition band)
%                   flt_pipeline(X,...,'iir',{[0.5 1.5],'highpass'},...)
%
%                   Example: implement an 1 Hz high-pass filter, passing the flt_iir arguments by name
%                   flt_pipeline(X,...,'iir',{'Frequencies',[0.5 1.5], 'Mode','highpass'},...)
%
%
%                'FIRFilter'/'fir' (see flt_fir): 
%                   apply an FIR-based frequency filter
%
%                   Example: implement an 1 Hz minimum-phase high-pass filter
%                   flt_pipeline(X,...,'fir',{[0.5 1.5],'highpass','minimum-phase'}, ...)
%
%                'Standardization'/'standardize' (see flt_standardize): 
%                   standardize the channels using a window of past signal
%
%                   Example: standardize using a moving 60-second window
%                   flt_pipeline(X,...,'Standardization',60,...)
%
%                'EpochExtraction'/'makepos'/'epoch' (see set_makepos): 
%                   Extract epochs from the signal (and deduce a target variable from the time-
%                   locking event types) (see set_makepos).
%
%                   Example: extract epochs around each occurrence of marker 'keypress1' and 
%                            'keypress2', cutting out segments that begin 2s before the marker and
%                            end 1s after the marker (note: keypress1 epochs will be assigned class 1
%                            and keypress2 epochs are being assigned class 2)
%                   flt_pipeline(X,...,'epoch',{[-2 1],{'keypress1','keypress2'}},...)
%
%                   Example: as before, but assign class 1 to keypress2 and class 2 to keypress1
%                   flt_pipeline(X,...,'epoch',{[-2 1],{'keypress2','keypress1'}},...)
%
%                   Example: assign class 1 to epochs around marker 'stimulus', and class 2 to epochs
%                            around markers 'keypress1' or 'keypress2'
%                   flt_pipeline(X,...,'epoch',{[-2 1],{'stimulus',{'keypress2','keypress1'}}},...)
%           
%                'BaselineRemoval'/'rmbase'/'baseline' (see flt_rmbase): 
%                   remove a baseline window of the given epoch(s)
%
%                   Example: for each epoch, average the signal value within -250ms to +100ms around 
%                            the time-locking event, and subtract that value from the entire epoch
%                            for the respective channel
%                   flt_pipeline(X,...,'BaselineRemoval',[-0.25 0.1],...)
%
%                'WindowSelection'/'window' (see flt_window): 
%                   apply a window function
%                   
%                   Example: scale each epoch in the data by a hann time window 
%                   flt_pipeline(X,...,'window','hann',...), 
%
%                   Example: restrict each epoch to the interval within 0 to 0.5 seconds after the 
%                            respective time-locking event
%                   flt_pipeline(X,...,'window',[0 0.5],...)
%
%                'SpectralSelection'/'spectrum' (see flt_spectrum): 
%                   apply a free-form spectral filter per epoch (note: if no epochs are present, this 
%                   gives a non-causal filtering of the entire signal)
%
%                   Example: apply a 7-30 Hz band-pass filter with linear falloffs at both edges
%                   flt_pipeline(X,...,'spectrum',[6.5 7.5 27 33],...)
%
%                'SparseReconstruction'/'reconstruct' (see flt_reconstruct): 
%                   reconstruct the signal in terms of a new (possibly overcomplete) basis
%
%                   Example: reconstruct the data in a random overcomplete basis using a fast EM method
%                   flt_pipeline(X,...,'reconstruct',{randn(32,1000), 'variant','FastEM'},...)
%
%                'ProjectionMatrix'/'project' (see flt_project): 
%                   apply a custom spatial projection matrix to the signal
%
%                   Example: project the data according to a random matrix
%                   flt_pipeline(X,...,'project',randn(128),...)
%
%                'SpectralTransform'/'fourier' (see flt_fourier):
%                   transform (usually epoched) signal into a Fourier representation 
%
%                   Example: flt_pipeline(X,...,'fourier','amplitude',...)
%
%
%                In addition, flt_pipeline provides afew of its own special-purpose parameters which
%                are not pipeline stages by themselves:
%
%                'FilterOrdering' : cell string array to partially override the default order
%                                   for the given pipeline stage functions. (e.g. {'set_makepos','flt_iir'})
%
%
% Out:
%   Signal  : the processed data set
%
% Examples:
%   % resample, apply an IIR high-pass filter and extract epochs (here: using short param names)
%   eeg = flt_pipeline(eeg, 'resample',200, 'iir',{[0.5 1],'highpass'}, 'epoch',[-1 1]);
%
%   % apply a surface laplacian and do an ICA decomposition, keeping the defaults for both stages
%   % (here: using the long parameter names)
%   eeg = flt_pipeline(eeg, 'SurfaceLaplacian',{}, 'ICA',{});
%
%   % as before, but use the 'on' syntax, instead of the cell array syntax
%   eeg = flt_pipeline(eeg, 'SurfaceLaplacian','on', 'ICA','on');
%
%   % assuming a function that runs flt_pipeline, but pre-specifies some of its own defaults internally,
%   % while allowing to customize/override all of its pipeline stage defaults, override the IIR filter parameters,
%   % and turn off the ICA decomposition (assuming that it was on as per defaults of specialpipeline())
%   eeg = specialpipeline(eeg, 'IIRFilter',{'[0.5 2],'highpass'}, 'ICA','off')
%
%   % as before, but use the [] syntax to turn off a pipeline stage
%   eeg = specialpipeline(eeg, 'IIRFilter',{'[0.5 2],'highpass'}, 'ICA',[])
%
%   % note: specialpipeline may be implemented as follows:
%   function signal = specialpipeline(signal,varargin)
%   ...
%   signal = flt_pipeline(signal,'IIRFilter',[5 7 25 30], 'ICA','on', varargin{:})
%   ...
%
%   % special: re-parse the list of supported plugin functions (affects GUIs displayed for flt_pipeline)
%   flt_pipeline('update')
%
% See also:
%   flt_clean_channels, flt_clean_peaks, flt_clean_spikes, flt_clean_windows (built-in Artifact Rejection)
%   flt_bandpower, flt_coherence, flt_fft, flt_fourier (built-in Spectral Transforms)
%   flt_fir, flt_iir, flt_spectrum (built-in Spectral Filtering)
%   flt_eog, flt_ica, flt_laplace, flt_project, flt_reref, flt_selchans, flt_seltypes, flt_selvolume, flt_stationary (built-in Spatial Filtering)
%   flt_epochica, flt_epochpca, flt_rmbase, flt_wavelet, flt_window (built-in Temporal Filtering)
%   flt_resample, flt_standardize, set_makepos, set_fit_dipoles (built-in Miscellaneous Filtering)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-29
dp;

warning off BCILAB:unresolved_dep

if ~isequal(varargin,{'update'})
    
    % --- collect, order and parameterize filters ---
    
    % scan the directories for available filters and assemble them in a list
    allflt = hlp_scope({'disable_dbstop_if_error_msg',true},@list_filters);
    
    % sort filters in the preferred order (and remove filters for which no preferences are known)
    custom_order = arg_extract(varargin,{'fltorder','FilterOrdering','FilterOrder','order'},[],{});
    [ordering,unlinked] = hlp_microcache('ordering',@order_filters,struct('name',{allflt.name},'properties',{allflt.properties}),custom_order);
    filters = allflt(ordering);
    
    % define/expose arguments to flt_pipeline (and define an argument for each filter)
    args = arg_define(@check_arguments, varargin, ...
        arg_norep({'signal','Signal'}), ...
        arg({'fltorder','FilterOrdering','order'},{},[],'Override filter order. Filters listed in this cell-string array are (partially) reordered according to this list. Example: {''set_makepos'',''flt_ica'',''flt_resample''}. See also the help of flt_pipeline.','type','cellstr','shape','row'), ...
        ... % list the argument specifications for all filters
        filters.spec);
    
    if any(cellfun(@(f)any(f==' '),args.fltorder))
        error('The FilterOrdering parameter is malformed; should be a cell array of strings, but was: %s',hlp_tostring(args.fltorder,10000)); end
    
    % check if there were conflicts / problems
    if ~isempty(unlinked)
        % there were ordering conflicts which led to the removal of filters from the chain this
        % should be extremely rare, if it can happen at all
        problematic = [];
        % check if some of those filters are actually in use during this invocation of the filter
        % pipeline
        for u=unlinked
            if args.(allflt(u).tag).arg_selection
                problematic(end+1) = u; end
        end
        if ~isempty(problematic)
            % if so, generate an error
            error('BCILAB:flt_pipeline:undefined_order',['The ordering relationship ' hlp_tostring(args.fltorder) ' is in conflict with the preferred ordering for the following filter stages:\n', ...
                hlp_tostring({allflt(probleatic).name}) '\n' ...
                'Their position in the filter chain is therefore not defined, but they are in use. Please include these filters in the filtering order definition, to specify their position.']);
        end
        % also note other conflicts (for filters that were not used in the current pipeline)
        disp(['Note: the ordering relationship ' hlp_tostring(args.fltorder) ' is in conflict with the ordering preferences for the following' quickif(~isempty(problematic),' further ',' ') 'nodes:']);
        disp(hlp_tostring({allflt(setdiff(unlinked,problematic)).name}));
        disp('Their position is currently not defined, but it is advisable to define it, if the filters are to be enabled eventually.');
    end

    
    % --- apply filters, as sig = flt_***(args,'signal',sig); ---
    
    for f=1:length(filters)
        % some new filters may have been added since the args were generated (thus we need to check)
        if isfield(args,['p' filters(f).tag])            
            setup = args.(['p' filters(f).tag]);
            if setup.arg_selection
                % also check if its dependencies are resolved
                for d = filters(f).properties.depends
                    if ~args.(['p' d{1}(5:end)]).arg_selection
                        warning('BCILAB:unresolved_dep',['Filter ' filters(f).name ' depends on ' d{1} ', which is currently disabled.']); end
                end
                args.signal = feval(filters(f).name,setup,'signal',args.signal);
            end
        end
    end

    signal = args.signal;
else
    % just update the list of filters and return
    hlp_scope({'disable_dbstop_if_error_msg',true},@list_filters,true);
end



% get a list of all existing filter specifications
function filters = list_filters(update_list)
dp;
debug = false;
% determine list of filters that are known to be unlistable here
known_incompliant = {'set_gettarget','set_combine','set_merge','set_joinepos','set_concat'};
% if SIFT is not installed, some more filters fall into that category
if ~exist('hlp_getModelingApproaches','file')
    known_incompliant = [known_incompliant {'flt_siftpipeline','flt_zscore'}]; end

persistent memo;
% if we need to (re-)collect the list
if isempty(memo) || exist('update_list','var') && update_list
    % get all directory entries
    found = [dir(env_translatepath('functions:/filters/flt_*.m')); 
             dir(env_translatepath('functions:/dataset_editing/set_*.m'));
             dir(env_translatepath('functions:/filters/in_development/flt_*.m'));
             dir(env_translatepath('functions:/dataset_editing/in_development/set_*.m'));
             dir(env_translatepath('private:/code/filters/flt_*.m'));
             dir(env_translatepath('private:/code/dataset_editing/set_*.m'));
             dir(env_translatepath('private:/code/filters/in_development/flt_*.m'));
             dir(env_translatepath('private:/code/dataset_editing/in_development/set_*.m'));
             dir(env_translatepath('home:/.bcilab/code/filters/flt_*.m'));
             dir(env_translatepath('home:/.bcilab/code/dataset_editing/set_*.m'))];

    % get file names, function names, function handles, function tags (names without prefix)
    files = setdiff({found.name},{'flt_pipeline.m','set_new.m','set_chanid.m','set_infer_chanlocs.m'});
    names = cellfun(@(n) n(1:end-2),files,'UniformOutput',false);
    funcs = cellfun(@(n) str2func(n),names,'UniformOutput',false);
    tags = cellfun(@(n) n(5:end),names,'UniformOutput',false);
    % obtain declared filter properties
    retain = true(size(funcs));
    for f=1:length(funcs)
        try
            props{f} = hlp_microcache('props',@arg_report,'properties',funcs{f},{hlp_fileinfo(funcs{f},[],'hash')});
        catch e
            if disp_once(['Cannot query properties of filter ' char(funcs{f}) ': ' e.message])  && debug
                hlp_handleerror(e); end
            retain(f) = false;
            continue;
        end
        % check/add missing properties
        for n = {'name','depends','precedes','follows','cannot_precede','cannot_follow'}
            prop = n{1};
            if ~isfield(props{f},prop)
                props{f}.(prop) = {}; end
            if ischar(props{f}.(prop))
                props{f}.(prop) = {props{f}.(prop)}; end
            if ~iscellstr(props{f}.(prop))
                error(['The given field must be a string or cell-string array: ' prop]); end
        end
        if ~isfield(props{f},'deprecated')
            props{f}.deprecated = false; end
        if ~isfield(props{f},'experimental')
            props{f}.experimental = false; end
        if ~isfield(props{f},'guru')
            props{f}.guru = false; end
    end
    
    % create argument specifications...
    for f=find(retain)
        % take the tag as the argument's code name, and the 'name' property as additional (e.g. GUI)
        % names
        nameset = [{['p' tags{f}]} props{f}.name tags(f)];
        % take the first line of the function's help text as description of the argument, followed
        % by the function name
        description = ['Implemented by ' names{f} '. ' hlp_fileinfo(names{f},[],'H1Line')];
        try
            specs{f} = arg_subtoggle(nameset,[],funcs{f},description,'deprecated',props{f}.deprecated,'experimental',props{f}.experimental); %#ok<*AGROW>
        catch e
            disp_once(['Cannot define argument slot for function ' char(funcs{f}) ' (likely an issue with the properties declaration): ' e.message]);
            retain(f) = false;
            continue;
        end
        if ~any(strcmp(char(funcs{f}),known_incompliant))
            if strcmp(debug,'hard')
                report = arg_report('rich',funcs{f}); %#ok<NASGU>
            else
                try
                    report = arg_report('rich',funcs{f}); %#ok<NASGU>
                catch e
                    % otherwise there is an actual error            
                    if ~disp_once(['Cannot query arguments of function ' char(funcs{f}) ' (likely an issue with the argument definition): ' e.message]) && debug
                        hlp_handleerror(e); end
                    retain(f) = false;
                    continue;
                end
            end
        else
            retain(f) = false;
        end
    end
    memo = struct('name',names(retain),'tag',tags(retain),'func',funcs(retain),'properties',props(retain),'spec',specs(retain));
end
% remember the result for the next time...
filters = memo;




% check the arguments passed into the function and optionally drop those that are not 
% known to flt_pipeline
function [args,fmt] = check_arguments(args,spec)
if ~isempty(args)
    if isstruct(args{1}) && (all(isfield(args{1},{'head','parts'})) || all(isfield(args{1},{'data','srate','tracking'})))
        % a data set was passed in as first argument: therefore, there is 1 positional argument and the rest are NVPs
        positionals = args(1);
        nvps = args(2:end);
    else
        % all arguments must be name/value pairs
        positionals = {};
        nvps = args;
    end
    
    % flatten the NVPs
    nvps = flatten_structs(nvps);
    
    % get all admissible names
    known_identifiers = {spec.names};
    known_identifiers_flat = [known_identifiers{:}];
    
    % collect all NVPs that are not supported
    dropped_empty = [];
    dropped_nonempty = [];
    for k=1:2:length(nvps)
        if strcmp(nvps{k},'arg_direct')
            continue; end
        if ~any(strcmp(nvps{k},known_identifiers_flat))
            if isempty(nvps{k+1}) || isequal(nvps{k+1},'off') || (isfield(nvps{k+1},'arg_selection') && isscalar(nvps{k+1}) && nvps{k+1}.arg_selection==0)
                dropped_empty(end+1) = k;
            else
                dropped_nonempty(end+1) = k;
            end
        end
    end
    % inform the user about what will be removed
    if ~isempty(dropped_empty)
        disp_once(['flt_pipeline: Unknown filter arguments were passed in as empty: ' hlp_tostring(nvps(dropped_empty)) '; these are being ignored.']); end
    if ~isempty(dropped_nonempty)
        warn_once(['flt_pipeline: Unknown filter arguments were passed in with non-trivial settings and are being ignored: ' hlp_tostring(nvps(dropped_nonempty)) '.']); end

    % ... and remove them from the arguments
    dropped_all = [dropped_empty dropped_nonempty];    
    nvps([dropped_all dropped_all+1]) = [];
    
    % construct outputs
    fmt = length(positionals);
    args = [positionals nvps];
else
    fmt = [0 1];
end


% substitute any structs in place of a name-value pair into the name-value list
function args = flatten_structs(args)
k = 1;
while k <= length(args)
    if isstruct(args{k})
        tmp = [fieldnames(args{k}) struct2cell(args{k})]';
        args = [args(1:k-1) tmp(:)' args(k+1:end)];
        k = k+numel(tmp);
    else
        k = k+2;
    end
end




% determine the preferred order of the filters
function [order,unlinked_filters] = order_filters(filters,override)
% optionally, a part of the order can be overridden with a cell array of names (or a cell array of
% edges)
if ~exist('override','var')
    override = {}; end

% we consider only filters for which we have some ordering information ('depends', 'follows', or
% 'precedes')
order_names = {};
for f=1:length(filters)
    flt = filters(f);
    % get a list of all filters related to flt
    related_filters = [flt.properties.depends flt.properties.follows flt.properties.precedes flt.properties.cannot_precede flt.properties.cannot_follow];
    % add to the full list
    order_names = [order_names related_filters];
    % also add the filter itself, if it relates to some other filter(s)
    if ~isempty(related_filters)
        order_names{end+1} = flt.name; end
end
order_names = unique(order_names);

if ~isempty(override)
    % also append those names that were specified in the override listing
    if iscellstr(override)
        order_names = [order_names override];
    elseif all(cellfun('isclass',override,'cell'))
        order_names = [order_names override{:}];
    else
        error('The override order of filters must either be a cell-string array {flt_xxx,flt_yyy,flt_zzz, ...} or a cell array of 2-element cell-string arrays {{flt_xxx,flt_yyy},{flt_xxx,flt_zzz}}');
    end
end

% retain only those filters that apppear in order_names
retain = [];
for f=1:length(filters)
    if any(strcmp(filters(f).name,order_names))
        retain(end+1) = f; end
end
filters = filters(retain);

% create a mapping from name to index in filters...
remap = struct();
for f=1:length(filters)
    remap.(filters(f).name) = f;  end
for f=1:length(filters)
    if isfield(filters(f),'properties') && isfield(filters(f).properties,'name')
        for n=1:length(filters(f).properties.name)
            name = filters(f).properties.name{n};
            if isfield(remap,name)
                error('BCILAB:flt_pipeline:duplicate_name','The given filter %s declares a name property (%s) that was already declared by some other filter.',filters(f).name,name);  end
            remap.(name) = f; 
        end
    end
end

% create a graph according to the filters' ordering preferences/constraints, as an edge list (of the
% type 'src precedes dst')
preferences = {};
constraints = {};
for i = 1:length(filters)
    flt = filters(i);
    for j = 1:length(flt.properties.depends)
        constraints{end+1} = [remap.(flt.properties.depends{j}),i]; end
    for j = 1:length(flt.properties.follows)
        preferences{end+1} = [remap.(flt.properties.follows{j}),i]; end
    for j = 1:length(flt.properties.precedes)
        preferences{end+1} = [i,remap.(flt.properties.precedes{j})]; end
    for j = 1:length(flt.properties.cannot_follow)
        constraints{end+1} = [i,remap.(flt.properties.cannot_follow{j})]; end
    for j = 1:length(flt.properties.cannot_precede)
        constraints{end+1} = [remap.(flt.properties.cannot_precede{j}),i]; end
end
edges = [preferences constraints];

newedges = {};
if ~isempty(override)
    % form the transitive closure of the graph, using Warshall's algorithm
    % (http://www.cs.nmsu.edu/~ipivkina/TransClosure/index.html)
    A = zeros(length(filters));
    for e=edges
        edge = e{1};
        A(edge(1),edge(2)) = 1; 
    end
    for i = 1:length(filters)
        for j = 1:length(filters)
            t = A(i,j) == 1 & A(j,:) == 1;
            A(i,t) = 1;
        end
    end
    % convert back to edge list...
    edges = {};
    for i = 1:length(filters)
        for j = 1:length(filters)
            if A(i,j)
               edges{end+1} = [i,j]; end
        end
    end
    
    % add override relationships
    if iscellstr(override)
        % a total order over a subset of nodes was specified
        for i=1:length(override)
            for j=i+1:length(override)
                newedges{end+1} = [remap.(override{i}) remap.(override{j})]; end
        end
    elseif all(cellfun('isclass',override,'cell'))
        % a partial order over a subset of nodes was specified
        for i=1:length(override)
            newedges{end+1} = [remap.(override{i}{1}) remap.(override{i}{2})]; end
    end
    edges = [edges newedges];
end

% remove duplicate edges
[dummy,edgeretain] = unique(cellfun(@num2str,edges,'UniformOutput',false)); %#ok<ASGLU>
edges = edges(edgeretain);

if ~isempty(override)
    % find all strongly connected components (cycles) of the graph
    % these are the conflicting ordering preferences
    nodes = repmat(struct('in',{[]},'out',{[]},'index',{[]},'lowlink',{[]}),1,length(filters));
    for e=edges
        edge = e{1};
        nodes(edge(1)).out(end+1) = edge(2);
        nodes(edge(2)).in(end+1) = edge(1);
    end       
    % use Tarjan's algorithm for this
    % (http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm)
    adj = sparse(cellfun(@(x)x(1),edges),cellfun(@(x)x(2),edges),true,length(retain),length(retain));
    comps = tarjan(full(adj)); %#ok<ACCUM>
    % remove all conflicting preference relationships (graph edges) ...
    offending_edges = {};
    for e=edges
        edge = e{1};
        if comps(edge(1)) == comps(edge(2))
            offending_edges{end+1} = edge; end
    end
    % ... except for the overridden ones, and the dependency edges
    [dummy,idx] = setdiff(cellfun(@num2str,offending_edges,'UniformOutput',false),cellfun(@num2str,[newedges constraints],'UniformOutput',false)); %#ok<ASGLU>
    edges_to_remove = offending_edges(idx);
    [dummy,edgeretain] = setdiff(cellfun(@num2str,edges,'UniformOutput',false),cellfun(@num2str,edges_to_remove,'UniformOutput',false)); %#ok<ASGLU>
    % note that this might leave some filter stages dangling in the air, as they might have been
    % part of some relationship that was in conflict with the override preferences if so, it is
    % strongly advisable to re-define their relationships...
    edges = edges(edgeretain);
    
    % re-check if we still have cycles in there...
    adj = sparse(cellfun(@(x)x(1),edges),cellfun(@(x)x(2),edges),true,length(retain),length(retain));
    compidx = tarjan(full(adj)); %#ok<ACCUM>
    node_count = accumarray(compidx',1,[max(compidx),1]);
    if any(node_count > 1)
        % yes, there are strongly connected components with more than 1 node
        problematic_components = {};
        for comp = find(node_count>1)
            problematic_components{end+1} = {filters(compidx == comp).name}; end
        % list them in an error message
        error('BCILAB:flt_pipeline:unresolvable',['The specified ordering is not compatible with the dependencies/constraints of the re-ordered filters.\nProblematic groups: ' hlp_tostring(problematic_components)]);
    end
end

% compute a topological order (http://en.wikipedia.org/wiki/Topological_sorting)
nodes = repmat(struct('in',{[]},'out',{[]}),1,length(retain));
for e=edges
    edge = e{1};
    nodes(edge(1)).out(end+1) = edge(2);
    nodes(edge(2)).in(end+1) = edge(1);
end

if ~isempty(override)
    % ... and also check whether the conflict resolution has erase all ordering relationships for 
    % some of the filters; we will warn about this in flt_pipeline
	unlinked_filters = find(cellfun('isempty',{nodes.in}) & cellfun('isempty',{nodes.out}));
    if ~isempty(unlinked_filters)
        % remove these from the ordered set
        idx_remap = 1:length(retain);
        idx_remap(unlinked_filters) = 0;
        idx_remap(find(idx_remap)) = 1:length(find(idx_remap)); %#ok<FNDSB>
        % remap the edge indices
        for e=1:length(edges)
            edges{e} = idx_remap(edges{e}); end
        retain(unlinked_filters) = [];
        % and recompute the node set
        nodes = repmat(struct('in',{[]},'out',{[]}),1,length(retain));
        for e=edges
            edge = e{1};
            nodes(edge(1)).out(end+1) = edge(2);
            nodes(edge(2)).in(end+1) = edge(1);
        end
    end
else
    unlinked_filters = [];
end

% create an index set over nodes with no in-edges
candidates = find(cellfun('isempty',{nodes.in}));
order = [];
while ~isempty(candidates)
    % take the first candidate
    [n,candidates] = deal(candidates(1),candidates(2:end));
    % put it into order
    order(end+1) = n;
    % for each node m with an edge from n to m...
    for m=nodes(n).out
        % remove the edge from the graph
        nodes(m).in = nodes(m).in(n ~= nodes(m).in);
        % if m has no other incoming edges, insert m into S
        if isempty(nodes(m).in)
            candidates(end+1) = m; end
    end
    nodes(n).out = [];
end
% return the ordered & reduced index set
order = retain(order);





