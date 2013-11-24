function should_evaluate = exp_beginfun(setting, varargin)
% Begin the definition of a function with lazy evaluation behavior.
% Evaluate = exp_beginfun(AttributeSetting, Attributes...)
%
% General
% =======
% If a function uses exp_beginfun/exp_endfun, MATLAB code that uses these functions does not "call"
% them in an imperative manner, but instead builds symbolic expressions that can be stored,
% transmitted, or manipulated before/if they are being evaluated. Evaluation is generally done via
% exp_eval() (or exp_eval_optimized()). Almost all uses of exp_beginfun by functions look like this:
%
%   function [res1,res2,...] = my_function(arg1,arg2, ...)
%   if ~exp_beginfun('somesetting') return; end
%
%   ... <actual implementation> ...
%
%   exp_endfun;
% 
% ... where 'somesetting' specifies the behavioral setting for the function. Generally, all
% functions that operate on data sets and modify the signal content in a way that needs to be
% replicated online have the setting 'filter', and all functions that have no equivalent during
% online processing (e.g., marker transcription) have the setting 'editing'. Functions that do
% processing that cannot be reproduced online have the setting 'offline'. See filters/flt_* and
% dataset_ops/set_* for some examples. A setting is a shortcut name for an assignment of attributes
% associated with the function declaration (attributes can also be directly specified to
% exp_beginfun).
%
%
% Filter and Editing Functions
% ============================
% Functions of type 'filter' should be written just like regular EEGLAB signal processing functions
% (i.e., they receive a data set and some options, and produce a transformed version of the data
% set), but functions that are supposed to be usable online should be "chunkable", i.e. they should
% produce the same sequence of output samples when they are called for consecutive segments (chunks)
% of the data set. Thus, the function should dump its internal state in an additional output (the
% 'state' variable) and be able to resume from this state when called on the next chunk. Chunks may
% well be just a single sample long. During online processing, state will be passed in automatically.
%
% In some cases, this strict behavior would be prohibitively inefficient - in this case, the
% function may refrain to non-adaptive behavior if the chunks are too short. In particular, in
% BCILAB, it may be assumed that each filter is first called on a complete calibration data set in a
% single call, and is then either called on a complete test data set in one call (the next 'chunk'),
% or piecemeal on small blocks of an incoming online stream. ICA, for example, adapts itself only on
% the calibration data set (i.e. the first chunk) and uses that decomposition for all remaining ones
% (i.e. during online use).
%
% Functions of type 'editing' should be written just like regular EEGLAB data set editing functions
% (e.g. pop_chanedit), BUT they must not change the data or srate fields.
%
% Rationale/Background
% ====================
% For the purpose of data processing function declarations, the exp_beginfun/exp_endfun contract
% handles the book-keeping necessary to be able to reproduce the filter (and entire filter chain)
% 1:1 online, and offers additional functionality, such as transparent caching of results (for later
% near-instantaneuous retrieval) and the annotation of all data sets with an executable expression
% which, when invoked, reproduces the data set 1:1, as well as checks to guarantee that no invalid
% "tampering" is done by the user on data sets (such as running a forward-backward filter in EEGLAB 
% or manually scaling the data). 
%
% In its general form, the functions exp_beginfun/exp_endfun are part of the machinery that enables
% symbolic computation in MATLAB (in the style of Mathematica), and the
% 'filter'/'editing'/'symbolic'/'annotate'/'default' argument is a short-cut which selects default
% attributes of the function being defined (which can be overridden). Some of these attributes are
% customizable plugin "slots" (plugins to be applied to each argument or all inputs or all outputs
% of the function being defined). 'Filter'/'editing' imply such hooks as default attributes, which
% in turn implement the tamper-proofing checks and some post-processing on the emitted data sets.
% The post-processing adds fields (subfields of the .tracking field) to the (otherwise arbitrary)
% data set structs that are produced by the function being defined, most importantly the
% .tracking.expression, which is an executable (recursive) expression representation which mirrors
% each data set in symbolic form, and .tracking.online_expression, in the case of filters, which
% defines the arguments to be used for the respective filter during online processing. These
% operations are usually completely transparent to the user (and locked away in the .tracking
% field), and make it exceptionally simple to author components for the toolbox.
%
% The symbolic-computation functionality implemented by the expressions/exp_* functions is largely
% non-overlapping with what the MATLAB Symbolic toolbox offers: the Symbolic toolbox exclusively
% implements mathematics (e.g. derivatives, integrals, simplifications), while the expression system
% of BCILAB exclusively implements expression manipulation, such as substitution, pattern matching,
% hashing, memoization, and reflection. These facilities allow to implement very high-level behavior
% in functions such as bci_train, such as generic programming (e.g. utl_crossval, utl_searchmodel,
% utl_gridsearch being non-intrusively adapted to EEG sets, STUDY sets, etc.), optimizations (e.g.
% in-memory caching and on-disk caching of intermediate results, computational shortcuts, and more
% in the future), with relatively little coding effort.
%
% While, by default, a function is invoked online with the same parameters that were used for
% offline processing (aside from a 'state' parameter for stateful functions), it is is possible to
% specify that certain parameters should have different values during online processing. This is
% enabled by the attribute 'append_online' (list of arguments to be appended in the online case),
% which can be specified either to exp_beginfun or exp_endfun. Alternatively, the entire expression
% can be subsituted by setting the 'set_online' attribute (expression to use during online
% processing).
%
% In:
%  AttributeSetting: baseline attribute assignment for the function; must be specified. Below are
%                    typical scenarios:
%                    * defining a function that should not evaluate immediately, but only after
%                      evaluation by an exp_eval() clause:
%                      'default' : delayed execution, input fingerprint checks
%
%                    * defining a function that processes data sets and should have
%                      delayed-evaluation semantics:
%                      'filter'  : functions that do non-trivial processing of data sets, and which
%                                  should support online replication of some offline processing;
%                                  attributes: delayed (1), add_impure (1), set_online
%                                  ('reproduce'), also performs EEGLAB-aware input data set checks &
%                                  transforms note: when an offline processing is replicated online,
%                                  the function is called successively on short EEGLAB data sets,
%                                  with inputs identical to the input parameters as they were
%                                  snapshot after execution of the function's body (by exp_endfun)
%                      'editing' : functions that do processing of data sets which has no equivalent
%                                  in online processing, i.e. which pass the input data through
%                                  unmodified; attributes: delayed (1), add_impure (1), set_online
%                                  ('passthrough'), also performs EEGLAB-aware input data set checks
%                      'offline' : functions that do processing of data sets which cannot be
%                                  reproduced online; these generate an error, when used within
%                                  online BCIss; attributes: delayed (1), add_impure (1), set_online
%                                  ('inapplicable'), also performs EEGLAB-aware input data set
%                                  checks
%
%                    * defining a function that performs purely symbolic manipulations of
%                      expressions:
%                      'symbolic': all inputs are held unevaluated, associated (impure) values are
%                                  ignored (thus no fingerprinting and no memoization); attributes:
%                                  delayed (1), hold ('all'), fingerprints unchecked & passed
%                                  through, memoize (0)
%                                  note: symbolic functions cannot be used in online code.
%
%                    * defining a function that performs immediate (non-lazy) processing of data
%                      structures, and which shall have their results annotated with an executable
%                      (repeatable) description of the processing performed:
%                      'annotate': delayed (0), add_impure (1)
%
%
%   Attributes... : optional name-value pairs to refine attributes of the function being defined:
%
%                     'delayed'   : whether the function has delayed-evaluation semantics, i.e. it 
%                                   doesn't evaluate unless passed to exp_eval (default: 1)
%
%                     'hold'      : []/'all'/'first'/'rest'; whether to leave all/the first/the
%                                   remaining arguments or no arguments of the function unevaluated
%                                   before executing the function body (useful for higher-order
%                                   functions), (default: [])
%
%                     'symbolic'  : whether the function being defined should be evaluated even if it
%                                   receives unevaluated expressions as arguments (ignoring those
%                                   arguments that fall under hold) - this is for symbol-manipulating
%                                   functions that can process unevaluated expressions (0/1, default: 0)
%
%                     'add_impure': integrate the original expression into the function value (when
%                                   evaluated), creating an 'impure expression', 0/1
%
%                     'set_online': Relevant only when the settings attribute is 'filter', 'editing'
%                                   or 'offline', may also be specified in exp_endfun. Introduces an
%                                   'online expression' into the function value (when evaluated),
%                                   which is a record of how the resulting data (usually a signal)
%                                   should be calculated online. Usually, this is identical to how
%                                   the filter is applied offline.
%                                   * 'passthrough': the function is skipped in the online
%                                                    processing
%                                   * 'reproduce': the function is called online with the same arguments
%                                                  as in offline processing
%                                   * 'inapplicable': the function gives an error when used in
%                                                     online processing
%                                   * {arg1, arg2, arg3, ...}: the given arguments will be used as
%                                                              function arguments during online use
%                                   * expression: the given expression will be used during online
%                                                 reproduction of the filter
%
%                     'append_online' : {arg1, arg2, arg3, ...}: the given arguments will be
%                                       appended to the function arguments during online use
%                                       (default: [])
%
%                     'argsteps', 'presteps', 'poststeps' : advanced attributes, see also "Advanced 
%                                                           Functionality" section.
%
%                     'delayed_online' : whether this function has delayed-evaluation semantics 
%                                        during online processing (default: false)
%
% Out:
%   Evaluate     : whether the body of the function being defined shall be evaluated or skipped
%
%
% Advanced functionality:
%   The allowed settings are customizable/extensible via the function exp_settings, which maps a
%   setting string onto a cell array of default attributes (name-value pairs).
%
%   During evaluation, each input argument is sent though a customizable sequence of functions,
%   which can be specified via the 'argsteps' attribute (as a cell array). Further, the attriutes
%   'presteps' and 'poststeps' are cell arrays of functions that should be invoked before and after
%   the function body runs, respectively; they operate on a context structure, which is defined
%   further below. These built-in attributes are the primary mechanism by which custom exp_settings 
%   add functionality to exp_beginfun/exp_endfun (e.g., annotating a dataset with an expression or
%   checking fingerprints).
%
%   The function's behavior can be controlled via dynamically scoped variables, as set by exp_block,
%   exp_set and exp_setdelayed, or correspondent attributes in exp_beginfun.
%     These variables can be assigned values, such as 0/1, but may be arbitrary expressions
%      * fingerprint_check: toggle per-argument fingerprint checking (default: true) - may use
%        @expression to refer to the argument expression, to implement expression-dependent checking
%      * fingerprint_create: toggle fingerprint creation (default: true) - may use @expression to
%        refer to the unevaluated complete expression, to implement expression-dependent checking
%      * memoize: selectively enable memoization; this is a cell array of the form
%        {location,expression,location,expression, ...}. location can be either 'disk' or 'memory',
%        expression should evaluate to true/false, may use @expression to refer to the unevaluated
%        complete expression, to implement expression-dependent memoization.
%
% See also:
%   exp_endfun, exp_eval, exp_settings
%
% Examples:
%   % a typical function declaration using the expression system
%   function result = myfunc(myarg,myotherarg,yetanotherarg)
%   if ~exp_beginfun('somesetting'), return; end
%   ...
%   exp_endfun;
%
%
%   % define a filter function using the expression system
%   % (note that, for filters to work with all BCILAB facilities, the varargin arguments should also 
%   %  be processed via arg_define)
%   function signal = flt_myfilter(varargin)
%   if ~exp_beginfun('filter'), return; end
%   ...
%   exp_endfun;
%
%
%   % as before, but define that the filter, when applied during online processing, should be passed
%   % a particular list of additional arguments
%   function signal = flt_myfilter(varargin)
%   if ~exp_beginfun('filter'), return; end
%   ...
%   exp_endfun('add_online',{'myinternalarg',10, 'my_reserved_argument','test'})
%
%
%   % like before, but this time replace the entire set of arguments; the assumption here is that
%   % the variable called signal contains the input signal, and that it can be passed to flt_myfilter 
%   % by name
%   function signal = flt_myfilter(varargin)
%   if ~exp_beginfun('filter'), return; end
%   ...
%   exp_endfun('set_online',{'signal',signal, 'myparam',10,'myotherparam',200})
%
%
%   % like before, but this time define that the filter, when applied online, should go through an
%   % entirely different function (here assumed to be defined using exp_beginfun itself) with specific
%   % arguments
%   function signal = flt_myfilter(varargin)
%   if ~exp_beginfun('filter'), return; end
%   ...
%   exp_endfun('set_online',flt_myfilter_online('signal',signal, 'specialparam1',10))
%
%
%   % define a dataset editing function
%   function signal = set_myoperation(varargin)
%   if ~exp_beginfun('editing'), return; end
%   ...
%   exp_endfun;
%
%
%   % define a function that operates on expressions symbolically
%   function result = my_operation(varargin)
%   if ~exp_beginfun('symbolic'), return; end
%   ...
%   exp_endfun;
%
%
%   % define a signal processing function that can only be applied in offline processing
%   function signal = set_myoperation(varargin)
%   if ~exp_beginfun('offline'), return; end
%   ...
%   exp_endfun;
%
%       
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-14

% --- obtain & interpret calling context ---

% obtain execution context
ctx.stack = dbstack;

% early check of whether the expression system is disabled
if  ~strcmp(setting,'symbolic') && hlp_resolve('disable_expressions',false,ctx) && ~isequal(varargin,{'delayed_online',true})
    % and exit exp_beginfun, if so
    assignin('caller','exp_internal_context',struct('enabled',0, 'exec_ctx',ctx));
    % the function should evaluate
    should_evaluate = 1;
    return;
end

% if we're being called by exp_eval, the function should be immediately evaluated - otherwise it may
% evaluate in a delayed fashion (i.e., only the expression itself is returned)
should_evaluate = any(strncmp('exp_eval',{ctx.stack(3:min(4,length(ctx.stack))).name},8));

% identify the calling function (i.e., the one using exp_beginfun)
callname = ctx.stack(2).name;
callfunc = str2func(callname);
[codehash,inargs,outargs] = utl_fileinfo([],callname); % note: elided the file for performance reasons


% --- read the caller's function arguments (unfolding varargs along the way) ---
try
    parts = {};
    for i=1:length(inargs)
        val = evalin('caller',inargs{i});
        if ~strcmp(inargs{i},'varargin')
            val = {val}; end
        parts = [parts val]; %#ok<AGROW>
    end
catch, end


% --- parse the attribute settings ---
if isempty(varargin) && strcmp(setting,'symbolic')
    % symbolic operations with no settings: fast path
    if ~should_evaluate
        % -> delayed evaluation: compose & return a new expression from transformed
        % arguments
        assignin('caller',outargs{1},struct('head',{callfunc},'parts',{parts}));
        % also assign [] to all other output arguments
        for k=2:length(outargs)
            assignin('caller',outargs{k},[]); end
    else
        % -> immediate evaluation: opt out of exp_endfun() processing (nothing to do for
        % symbolic functions)
        assignin('caller','exp_internal_context',struct('enabled',{0}, 'exec_ctx',{ctx}));
    end
    return;
else
    % everything else: regular path
    
    % define general attributes
    attribs = {'delayed',1,'add_impure',1,'hold',[],'symbolic',0,'set_online',[],'append_online',[],'delayed_online',false,'argsteps',[],'presteps',[],'poststeps',[],'fingerprint_check',[],'fingerprint_create',[],'memoize',[]};
    
    % set defaults according to the chosen setting
    switch setting
        case 'symbolic'
            defaults = {'add_impure',0, 'symbolic',1, 'hold','all', 'fingerprint_check',0, 'fingerprint_create','passthrough', 'memoize',0};
        case 'default'
            defaults = {'add_impure',0, 'argsteps',{@utl_check_fingerprint}};
        case 'annotate'
            defaults = {'delayed',0};
        otherwise
            defaults = exp_settings(setting);
    end
    
    % concatenate general attribs, per-setting defaults and user overrides into one list of
    % name-value pairs (NVPs)
    nvps = [attribs defaults varargin];
    % retain only the last assignment for each name
    [s,indices] = sort(nvps(1:2:end));
    indices(strcmp(s((1:end-1)'),s((2:end)'))) = [];    
    % and turn them into an options struct
    opts = cell2struct(nvps(2*indices),nvps(2*indices-1),2);
end

% if delayed is false, it enforces immediate evaluation
if ~opts.delayed
    should_evaluate = true; end


if should_evaluate

    % apply per-argument transformation steps (and pass the full options)
    expression = struct('head',{callfunc},'parts',{parts},'codehash',{codehash});
    for step = opts.argsteps
        for p=1:length(parts)
            parts{p} = step{1}(parts{p},opts,ctx,expression); end
    end

    % --- check if the result of evaluating this expression is already in the cache ---
      
    % purify the expression    
    expression = utl_purify_expression(struct('head',{callfunc},'parts',{parts},'codehash',{codehash}));
    
    % quick check if memoization is disabled
    if ~isequal(opts.memoize,0)
        % do a cache lookup
        [action,result] = utl_memoize_lookup(expression,opts.memoize,ctx);
        % act according to outcome
        if strcmp(action,'return')
            % assign the retrieved data to the output argument(s)
            for k=1:length(result)
                assignin('caller',outargs{k},result{k}); end
            % and make sure that the function does not evaluate further
            should_evaluate = 0;
            return;
        elseif strcmp(action,'memoize')
            % remember to memoize the expression later, under the given id(s)
            memoize_id = result;
        else
            % do not memoize the expression
            memoize_id = [];
        end
    else
        memoize_id = [];
    end
    
    
    % --- not in cache: evaluate arguments  ---
    
    % evaluate the function's arguments (recursively), if they are not held
    if ~isempty(parts) && ~any(strcmp(opts.hold,{'all','first'}));
        parts{1} = exp_eval(parts{1},inf);
        if ~opts.symbolic && ~is_evaluated(parts{1})
            % one of the inputs could not be evaluated (e.g. using an undefined symbol): return the
            % expression
            assignin('caller',outargs{1},expression);
            % also assign [] to all other output arguments
            for k=2:length(outargs)
                assignin('caller',outargs{k},[]); end
            should_evaluate = false;
            return;
        end
    end    
    if ~any(strcmp(opts.hold,{'all','rest'}))
        for i=2:length(parts)
            parts{i} = exp_eval(parts{i},inf); %#ok<AGROW>
            if ~opts.symbolic && ~is_evaluated(parts{i})
                % one of the inputs could not be evaluated (e.g. using an undefined symbol): return
                % the expression
                assignin('caller',outargs{1},expression);
                % also assign [] to all other output arguments
                for k=2:length(outargs)
                    assignin('caller',outargs{k},[]); end
                should_evaluate = false;
                return;
            end
        end
    end
    
    % --- prepare function execution ---
    
    % begin measuring the evaluation time of the expression ...
    eval_time = tic;
    % ... and measure the size of the inputs (both used to determine the best cache location)
    stats = whos('parts');
    input_size = stats.bytes;
    
    % create a context struct to keep track of the inputs, outputs, and attributes of the function
    % throughout its evaluation
    context = struct( ...
        'inargs',{inargs}, 'outargs',{outargs}, ...                     % input/output argument names
        'ws_input_pre',{struct()}, 'ws_output_post',{struct()}, ...     % input/output function workspace before/after evaluation of the function's body
        'expression_preeval',{expression}, 'expression_posteval', struct('head',{callfunc},'parts',{parts}), ... % expression arguments (original and evaluated)
        'opts',{opts}, 'eval_time',{eval_time}, 'input_size',{input_size}, 'memoize_id',{memoize_id}, 'exec_ctx',{ctx}, 'enabled',{1});    % misc fields
    
    % re-derive the function's input & output workspace from the now evaluated parts (folding
    % varargs again)
    for i=1:length(parts)
        if strcmp(inargs{i},'varargin')
            context.ws_input_pre.(inargs{i}) = parts(i:end);
            break;
        else
            context.ws_input_pre.(inargs{i}) = parts{i};
        end
    end
    
    % execute the prestep(s) on the function workspace
    for step = opts.presteps'
        context = step{1}(context); end
    
    % store the context object...
    assignin('caller','exp_internal_context',context);

    % assign all input for the function body's execution
    for fn=fieldnames(context.ws_input_pre)'
        assignin('caller',fn{1},context.ws_input_pre.(fn{1})); end
    
    % (from here on, the function body will execute)
else
    
    % --- delayed evaluation: return the expression ---

    expression = struct('head',{callfunc},'parts',{parts},'codehash',{codehash});
    assignin('caller',outargs{1},expression);
    % also assign [] to all other output arguments
    for k=2:length(outargs)
        assignin('caller',outargs{k},[]); end
    
    % (the function body will return immediately, as the should_evaluate output is false)
end
