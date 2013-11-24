function exp_endfun(varargin)
% End the definition of an expression-aware function.
% exp_endfun(Arguments...)
%
% See exp_beginfun. Invokes plugin functions that shall run after the body of the expression-aware 
% function (here called poststeps) and handles memoization.
%
% In:
%   Arguments... : optional name-value pairs; some options allowed in exp_beginfun can also be
%                  overridden here: 'set_online', 'append_online', 'add_impure'
%
% Notes:
%   The exp_endfun must be called in every code path that exits the function being defined (except
%   for error conditions), i.e. a premature return must be preceded by and exp_endfun;
%
% See also:
%   exp_beginfun, exp_eval
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15

% retrieve the context structure
% If you get an "Undefined function or variable" error here then you forgot to call exp_beginfun
% at the beginning of the same function in which exp_endfun was called
context = evalin('caller','exp_internal_context');
if ~context.enabled
    return; end

% update the new function input & output arguments
for i=1:length(context.outargs)
    try context.ws_output_post.(context.outargs{i}) = evalin('caller',context.outargs{i}); catch, end
end

% if the arg system is used by the function, we have an easier-to-handle version of the posteval expression
if evalin('caller','exist(''arg_nvps'',''var'')')
    context.expression_posteval.parts = evalin('caller','arg_nvps'); end

% assign additional exp_endfun options to the context's options...
for k=1:2:length(varargin)
    context.opts.(varargin{k}) = varargin{k+1}; end

if context.opts.add_impure
    % if impure, incorporate the original unevaluated expression into the function result (if it's a struct)
    if isfield(context.ws_output_post,context.outargs{1}) && isstruct(context.ws_output_post.(context.outargs{1}))
        context.ws_output_post.(context.outargs{1}).tracking.expression = context.expression_preeval; end
end

% run the poststep(s) on the workspaces (add online expression, etc.)
for s = 1:length(context.opts.poststeps)
    context = context.opts.poststeps{s}(context); end

% handle fingerprinting of the result
if isempty(context.opts.fingerprint_create)
    context.opts.fingerprint_create = hlp_resolve('fingerprint_create',true,context.exec_ctx); end
if strcmp(context.opts.fingerprint_create,'passthrough')
    % no changes are being made to the fingerprint of the output
elseif isfield(context.ws_output_post,context.outargs{1}) && is_impure_expression(context.ws_output_post.(context.outargs{1}))
    % we have an impure expression as output
    if ~isequal(false,context.opts.fingerprint_create) && ~isequal(false,exp_eval(utl_replacerepeated(context.opts.fingerprint_create,{exp_rule(@expression,context.expression_preeval)}),inf)) ...
        % fingerprint creation turned on: fingerprint the data, so we can later check the consistency of impure expression
        context.ws_output_post.(context.outargs{1}).tracking.fingerprint = hlp_fingerprint(rmfield(context.ws_output_post.(context.outargs{1}),'tracking'));
    else
        % fingerprint creation turned off: make sure that the fingerprint gets removed...
        if isfield(context.ws_output_post.(context.outargs{1}).tracking,'fingerprint')
            context.ws_output_post.(context.outargs{1}).tracking = rmfield(context.ws_output_post.(context.outargs{1}).tracking,'fingerprint'); end
    end
end

% remember the computation time of the object
if isfield(context.ws_output_post,context.outargs{1}) && is_impure_expression(context.ws_output_post.(context.outargs{1}))
    context.ws_output_post.(context.outargs{1}).tracking.computation_time = toc(context.eval_time); end

% assign the changed values in the workspace
for fn=fieldnames(context.ws_output_post)'
    assignin('caller',fn{1},context.ws_output_post.(fn{1})); end

% memoize the results (for the original expression)
if ~isempty(context.memoize_id) && isfield(context.ws_output_post,context.outargs{1})
    obj = {context.ws_output_post.(context.outargs{1})};
    for k=2:length(context.outargs)
        % and also append the other output arguments as far as they are set
        try obj = [obj {context.ws_output_post.(context.outargs{k})}]; catch break; end
    end
    
    % for each location... (currently: memory or disk)
    for c=1:length(context.memoize_id)
        id = context.memoize_id{c};
        utl_memoize_commit(obj,id,context.input_size);
    end
end
