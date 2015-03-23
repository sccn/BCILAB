function context = utl_add_online(context)
% Internal. Post-step for exp_beginfun, adding an online expression to a dataset.
%
% This is a poststep (i.e. calculation that runs at the time of the exp_endfun call) to insert a
% .tracking.online_expression field into the data set that is returned by a filter stage; the
% online_expression is the one that will be evaluated online per chunk to realize this filter stage.
%
% This is the function that handles the 'set_online' and 'append_online' attributes that can be 
% set via exp_beginfun and overridden in exp_endfun.
%
% In:
%   Context : The evaluation context struct used by exp_beginfun/exp_endfun. The relevant fields are
%             * .opts.set_online (the value of the set_online attribute)
%
%               The value of this attribute is determined by the setting passed to exp_beginfun,
%               where the 'filter' setting sets it to 'reproduce', the 'editing' setting sets it to
%               'passthrough', and the 'offline' setting sets it to 'inapplicable'; some filters 
%               override it according to the permitted syntax explained in exp_beginfun, which is 
%               either a cell array of new arguments to use, or an expression to use.
%
%             * .opts.append_online (the value of the set_online attribute); this attribute is
%               by default not set, unless a filter sets it explicitly. Some filters use it to append
%               extra arguments during online processing.
%
% Out:
%   Context : The updated evaluation context. This function will only update 
%             the .tracking.online_expression field of the first output of the function that uses
%             exp_endfun, which is held in a field of the context for later assignment to the 
%             function's workspace by exp_endfun; the context field is:
%             context.ws_output_post.(context.outargs{1}).tracking.online_expression
%
% See also:
%   exp_beginfun, exp_endfun, utl_complete_model


% --- handle the 'set_online' attribute ---

% we construct the final_expression based on the value of this attribute
set_online = context.opts.set_online;
if ischar(set_online)
    switch set_online
        case 'reproduce'
            % take the original expression directly as the online expression (this is the default for
            % 'filter' expressions, as specified in exp_settings)
            final_expression = context.expression_posteval;            
        case 'passthrough'
            % the function shall be skipped online: we effectively set as this online expression one of our
            % input signal's online expressions (the data sets among the filter arrguments are here viewed
            % as the input 'signals' during online operations). This is the default for the 'editing' setting 
            % in exp_begindef.
            online_expressions = {};
            for k=1:length(context.expression_posteval.parts)
                exp = context.expression_posteval.parts{k};
                if is_impure_expression(exp) && isfield(exp.tracking,'online_expression')
                    online_expressions{end+1} = exp.tracking.online_expression; end %#ok<AGROW>
            end

            if isempty(online_expressions)
                % there must be at least one signal among the inputs
                final_expression = struct('head',@error,'parts',{{'BCILAB:exp_beginfun:no_skip','This stage cannot be skipped online, since it does have not any input signal.'}});
            elseif length(online_expressions) > 1 && ~all(cellfun(@(e)utl_same(e,online_expressions{1}), online_expressions(2:end)))
                % and if there are multiple input signals, they must all be the same
                final_expression = struct('head',@error,'parts',{{'BCILAB:exp_beginfun:no_skip','This stage cannot be skipped online, since it has multiple different input signals.'}});
            else
                final_expression = online_expressions{1};
            end

            if ~isempty(context.opts.append_online)
                error('BCILAB:exp_beginfun:cannot_append','You cannot append expressions to a filter that has been declared as ''editing'' (i.e. skipped). However, you can replace the entire expression for it by what you like to invoke online using the ''set_online'' attribute.'); end
        case 'inapplicable'
            % Generate an error when used online. This is the default for the 'offline' setting in exp_begindef.
            final_expression = struct('head',@error,'parts',{{'BCILAB:exp_beginfun:no_online','This function cannot be run online.'}});
        case 'imported'
            % Set the initial online expression for data that was freshly imported
            res = context.ws_output_post.(context.outargs{1});            
            % check if the object already has an online expression (e.g., result of io_loadset or processed version thereof)
            if isfield(res,'tracking') && all(isfield(res.tracking,{'online_expression','fingerprint'}))
                if strcmp(char(raw.tracking.online_expression.head),'rawdata')
                    % if this is a rawdata (trivial) expression, we can do at least as good a job here if we
                    % override it based on the data that we actually got in res
                    final_expression = struct('head',@rawdata,'parts',{{{res.chanlocs.labels},unique({res.chanlocs.type})}});
                elseif isequal(hlp_fingerprint(rmfield(res,'tracking')),res.tracking.fingerprint)
                    % the data has been modified since the expression was evaluated: the online expression does not apply and we drop it
                    % this is not unusual, e.g., when the imported data was manually edited after an
                    % initial import and/or processing using io_loadset
                    final_expression = struct('head',@rawdata,'parts',{{{res.chanlocs.labels},unique({res.chanlocs.type})}});
                else
                    % the data has a non-trivial expression that matches: we keep it by notify the
                    % user that an importer is not the best place to do signal processing
                    disp_once('Best practice for an import function is to leave most signal processing to the BCI approach.');
                    final_expression = res.tracking.online_expression;
                end
            else
                % there is no online expression available, so we build it 
                final_expression = struct('head',@rawdata,'parts',{{{res.chanlocs.labels},unique({res.chanlocs.type})}});
            end
        otherwise
            error('Unsupported value for the ''set_online'' attribute of exp_beginfun/exp_endfun: %s',set_online);
    end    
elseif isfield(set_online,{'head','parts'})
    % completely replace the original expression by the set_online attribute
    final_expression = set_online;
elseif iscell(set_online)
    % leave the function the same and override only the arguments (parts) by by the set_online
    % attribute
    final_expression = context.expression_posteval;
    final_expression.parts = set_online;
elseif isa(set_online,'function_handle')
    % the given expression is just a function handle (i.e. a symbol)
    final_expression = set_online;
else
    % unknown format
    error('Unsupported format for the ''set_online'' attribute of exp_beginfun/exp_endfun: %s',hlp_tostring(set_online));
end

% --- handle the append_online option ---

% append additional parameters, if requested
if ~isempty(context.opts.append_online)
    if iscell(context.opts.append_online)
        final_expression.parts = [final_expression.parts context.opts.append_online];
    else
        error('The ''append_online'' attribute of exp_beginfun/exp_endfun must be set to a cell array of parameters, but is: %s',hlp_tostring(context.opts.append_online));
    end
end

% --- finalize the final_expression ---

% append the state to the expression
if length(context.outargs) > 1 && isfield(context.ws_output_post,context.outargs{2})
    % if a second output is present and assigned, we treat it as a state output, and include it
    % in the final expression this will be picked up by the online system
    final_expression.state = context.ws_output_post.(context.outargs{2});
    % also, generally append 'state',state to the expression - for proper cross-validation behavior
    final_expression.parts = [final_expression.parts {'state' final_expression.state}];
end
    
% for all sub-arguments that have an online expression (i.e. data sets / signals), substitute 
% their expression into the final online expression
if isfield(final_expression,'parts')
    for k=1:length(final_expression.parts)
        exp = final_expression.parts{k};
        if isfield(exp,'tracking') && isfield(exp.tracking,'online_expression')
            final_expression.parts{k} = exp.tracking.online_expression; end
    end
end
    
% finally append an arg_direct,1 to the final expression
final_expression.parts{end+1} = struct('arg_direct',{1});

% assign the expression to the .tracking.online_expression field of the first output (-to-be)
if isfield(context.ws_output_post.(context.outargs{1}),'tracking')
    context.ws_output_post.(context.outargs{1}).tracking.online_expression = final_expression; end

