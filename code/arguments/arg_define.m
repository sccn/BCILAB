function outstruct = arg_define(vals,varargin)
% Parse function arguments and enable function introspection.
% Struct = arg_define(Values, Specification...)
% Struct = arg_define(Format, Values, Specification...)
%
% Replacement for the parameter declaration line of a function. Typically that Function will only
% declare a single varargin argument and then calls arg_define(varargin, ...) to parse it into a
% Struct. The Specification of allowed arguments with defaults, range and help text is given as a
% as a series of arg(), arg_sub(), etc. definitions, each of which declares one named argument.
%
% One can also omit the output Struct to have the arguments assigned directly into the Function's
% workspace, but that is discouraged  as it will only work for argument names that do not clash with
% any function on the MATLAB path (due to a limitation in MATLAB).
%
% The format of varargin is flexible and can be a list of a fixed number of positional arguments
% (i.e., the typical MATLAB calling format), optionally followed by a list of name-value pairs
% (NVPs, e.g., as the format accepted by figure()). The name-value pairs may be interleaved with
% structs, i.e., one may pass a mix of 'name',value,STRUCT,'name',value,'name',value, ...), which is
% not ambiguous. The allowed number of positional arguments, if any, can be specified by the
% optional extra Format argument to arg_define. Only names that are listed in the Specification may
% be used as named arguments.
%
% A main feature of arg_define() is that the argument specification can be reported to outside
% functions (queried using arg_report()), which allows to auto-generate GUIs, help text, and various
% other things.
%
% In:
%   Format : Can be omitted. The number of allowed positional arguments of the Function, or a vector
%            of multiple possibilities. More precisely, this is the number of leading arguments that
%            can be passed to the Function as positional arguments, while the remaining arguments
%            are interpreted as a list of name-value pairs and/or structs. See Format Options below
%            for additional options and details. (default: [0 Inf])
%
%   Values : A cell array of values passed to the function (usually the calling function's
%            "varargin"). Interpreted according to the Format and the Specification.
%
%   Specification... : The specification of the calling function's arguments; this is a sequence of
%                      arg(), arg_norep(), arg_nogui(), arg_deprecated(), arg_sub(),
%                      arg_subswitch(), arg_subtoggle() specifiers.
%
% Out:
%   Struct : A struct with values assigned to fields, according to the Specification and Format.
%
%            Note: If this output is not requested by the Function, the contents of Struct are
%            instead assigned to the Function's workspace -- but note that this only works for
%            variable names are *not* also names of functions in the path (due to a deficiency in
%            MATLAB's treatment of variable identifiers). Thus, it is good advice to use variable
%            names that are unlikely to be function names to avoid this situation (e.g.,
%            long/expressive or CamelCase names).
%
% See also:
%   arg, arg_nogui, arg_norep, arg_deprecated, arg_sub, arg_subswitch, arg_subtoggle
%
% Format Options:
%    * If Format is a vector of multiple values, then the smallest number for which the 
%      Function's inputs are consistent, in the sense that each name in Values matches a
%      name in the Specification, is accepted.
%    * If Format is numeric but does not contain a 0, a 0 will be implicitly prepended
%      (that is, one can always pass in everything as name-value pairs or structs).
%    * If Format is [], any number of positional arguments is allowed (i.e., Format defaults
%      to 0:length(Specification)), although this is discouraged in practice as it may
%      result in mis-parsing when mis-spelled arguments are passed in (arg_define would then
%      assume that anything up to the mis-spelled argument is passed in positionally).
%    * If Format is a function handle, the given function can be used to transform the
%      Values prior to any other processing into a new Values cell array. The function may
%      optionally specify a new (numeric) Format as its second output argument (if not
%      specified, this is 0). The function may optionally take in the argument specification
%      (struct array of arg_specifier), although some details of this data structure are
%      internal and subject to change in future releases (particularly fields marked as
%      INTERNAL).
%
% Performance Tips:
%   1) If a name-value pair 'arg_direct',true is passed (preferably as last argument for best
%      performance), or a struct with a field named 'arg_direct' is passed (that is set to true),
%      then only the values that are passed in will be assigned (and all type checking and default
%      values are skipped); note that positional arguments are unsupported in this mode. This is a
%      fast way to call a function (not much slower than hand-written argument parsing), and it
%      usually means that all arguments should be passed in. To obtain a struct with the full set of
%      function arguments to pass in (and their defaults), one can use arg_report.
%
%      Please do not pass multiple occurrences of 'arg_direct' with conflicting values to
%      arg_define, since the resulting behavior is then undefined.
%
%   2) If a function has many arguments it is faster to receive the arguments in a struct (optional
%      return value of arg_define) rather than to have them assigned to the function workspace.
%
% Internal Arguments and Return Values:
%   The Function may be called with the request to deliver the parameter specification as opposed to
%   following the normal execution flow (this is done by arg_report), an effectively implemented by
%   passing in the extra name-value pair '__arg_report__',true at the end. arg_define responds to
%   this task by throwing an exception of the type 'BCILAB:arg:report_args' using arg_issuereport,
%   which is caught by arg_report to obtain the specification struct. When a report is issued the
%   function may be called with some further internal arguments preceding the '__arg_report__',true,
%   including __arg_skip__,true (skip leading skippable arguments) and __arg_nodefaults__,true (do
%   not apply default values).
%
% Examples:
%   function myfunction(varargin)
%
%   % begin a default argument declaration and declare a few arguments; The arguments can be passed either:
%   % - by position: myfunction(4,20); including the option to leave some values at their defaults, e.g. myfunction(4) or myfunction()
%   % - by name: myfunction('test',4,'blah',20); myfunction('blah',21,'test',4); myfunction('blah',22);
%   % - as a struct: myfunction(struct('test',4,'blah',20))
%   % - as a sequence of either name-value pairs or structs: myfunction('test',4,struct('blah',20)) (note that this is not ambiguous, as the struct would come in a place where only a name could show up otherwise
%   arg_define(varargin, ...
%       arg('test',3,[],'A test.'), ...
%       arg('blah',25,[],'Blah.'));
%
%   % a special syntax that is allowed is passing a particular parameter multiple times - in which case only the last specification is effective
%   % myfunction('test',11, 'blah',21, 'test',3, struct('blah',15,'test',5), 'test',10) --> test will be 10, blah will be 15
%
%   % begin an argument declaration which allows 0 positional arguments (i.e. everything must be passed by name
%   arg_define(0,varargin, ...
%
%   % begin an argument declaration which allows exactly 1 positional arguments, i.e. the first one must be passed by position and the other one by name (or struct)
%   % valid calls would be: myfunction(3,'blah',25); myfunction(3); myfunction(); (the last one assumes the default for both)
%   arg_define(1,varargin, ...
%       arg('test',3,[],'A test.'), ...
%       arg('blah',25,[],'Blah.'));
%
%   % begin an argument decalration which allows either 2 positional arguments or 0 positional arguments (i.e. either the first two are passed by position, or all are passed by name)
%   % some valid calls are: myfunction(4,20,'flag',true); myfunction(4,20); myfunction(4,20,'xyz','test','flag',true); myfunction(4); myfunction('flag',true,'test',4,'blah',21); myfunction('flag',true)
%   arg_define([0 2],varargin, ...
%       arg('test',3,[],'A test.'), ...
%       arg('blah',25,[],'Blah.'), ...
%       arg('xyz','defaultstr',[],'XYZ.'), ...
%       arg('flag',false,[],'Some flag.'));
%
%   % begin an argument declaration in which the formatting of arguments is completely arbitrary, and a custom function takes care of bringing them into a form understood by
%   % the arg_define implementation. This function takes a cell array of arguments (in any formatting), and returns a cell array of a standard formatting (e.g. name-value pairs, or structs)
%   arg_define(@myparser,varargin, ...
%       arg('test',3,[],'A test.'), ...
%       arg('blah',25,[],'Blah.'));
%
%   % return the arguments as fields in a struct (here: opts), instead of directly in the workspace
%   opts = arg_define(varargin, ...
%       arg('test',3,[],'A test.'), ...
%       arg('blah',25,[],'Blah.'));
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-24

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

    % first parse the inputs to arg_define
    [fmt,vals,compressed_spec,structmask,report_type,skip,nodefaults] = process_inputs(vals,varargin);
    
    % check if we are in direct mode (fast shortcut)
    if strcmp(report_type,'none') && is_direct_mode(vals,structmask)
        % in direct mode we assign only those variables that were passed in
        [nvps,outstruct] = assign_direct(vals,structmask,nargout);
    else
        % otherwise we perform full parsing
        
        % get the calling function's name (for cache lookups and diagnostics)
        caller_name = hlp_getcaller;
        if ~isvarname(caller_name)
            caller_name = 'anonymous'; end
        
        % expand specification into a struct array, derive some properties
        [spec,flat_names,first_names,name2idx,leading_skippable,checks] = process_spec_cached(caller_name,compressed_spec,report_type,~nodefaults,nargout==0);

        % convert vals to a canonical list of name-value pairs (NVPs)
        nvps = arguments_to_nvps(caller_name,fmt,vals,structmask,flat_names,first_names,skip*leading_skippable);

        % assign the NVPs to the spec
        spec = assign_nvps(spec,nvps,name2idx,report_type,caller_name,true);

        % optionally remove unassigned spec entries
        if nodefaults && ~isempty(spec)
            spec(strcmp('__arg_unassigned__',{spec.value})) = []; end
                
        % generate outputs from spec
        switch report_type
            case 'none'
                % reset the skippable flag (these arguments are not skippable as far as the Function 
                % to which we're returning them is concerned)
                [spec.skippable] = deal(false);
                % handle mandatory entries
                mandatory_entries = find(strcmp('__arg_mandatory__',{spec.value}));
                if mandatory_entries
                    error(['The arguments ' format_cellstr({spec(mandatory_entries).first_name}) ' were unspecified but are mandatory.']); end
                % build output struct and generate full NVP list for assignment to workspace
                outstruct = arg_tovals(spec,[],'struct',false,checks.unassigned,checks.expression,checks.conversion);
                nvps = reshape([fieldnames(outstruct)';struct2cell(outstruct)'],1,[]);
            case 'vals'
                arg_issuereport(arg_tovals(spec,[],'struct',false,checks.unassigned,checks.expression,checks.conversion));
            case 'nvps'
                arg_issuereport(arg_tovals(spec,[],'cell',false,checks.unassigned,checks.expression,checks.conversion));
            case {'lean','rich'}
                arg_issuereport(spec);
            otherwise
                error(['Unrecognized report type requested: ' report_type]);
        end
    end
    
    % if requested, place the arguments in the caller's workspace
    if ~nargout
        try
            for k=1:2:length(nvps)
                assignin('caller',nvps{k},nvps{k+1}); end
        catch e
            if strcmp(e.identifier,'MATLAB:err_static_workspace_violation')
                error('In a function with nested functions you need to capture the outputs of arg_define into a struct.');
            else
                rethrow(e);
            end
        end
    end

    % also generally return the arguments in their native NVP/struct form
    try        
        assignin('caller','arg_nvps',nvps);
    catch %#ok<CTCH>
        % this operation might be disallowed under some circumstances (e.g., arg_define in a nested function)
    end
end


% process the inputs to arg_define into Format, Values and Specification
% also precompute the structmask (bitmask that encodes which elements in Values are structs)
% and split both the report_type, skip and nodefaults flags off from vals
function [fmt,vals,compressed_spec,structmask,report_type,skip,nodefaults] = process_inputs(vals,compressed_spec)
    skip = false;
    nodefaults = false;
    if iscell(vals)
        % no Format specifier was given: use default
        fmt = [0 Inf];
    else
        % a Format specifier was given as first argument (need to shift remaining arguments by 1)
        fmt = vals;
        vals = compressed_spec{1};
        compressed_spec(1) = [];
    end
    
    % extract report type
    if length(vals)>1 && isequal(vals{end-1},'__arg_report__')
        report_type = vals{end};
        % check for some report types that can be handled immediately
        switch report_type
            case 'raw'
                arg_issuereport(compressed_spec);
            case 'properties'
                arg_issuereport(struct());
            case {'lean','rich'}
                if length(vals) == 2
                    arg_issuereport(hlp_microcache('spec',@expand_spec,compressed_spec,report_type,true,hlp_getcaller(2))); end
            case 'handle'
                error('To make function handles accessible, use the function expose_handles().');
            case 'supported'
                arg_issuereport(true);
        end
        vals(end-1:end) = [];
        
        % check for nodefaults flag (do not apply defaults)
        if length(vals)>1 && isequal(vals{end-1},'__arg_nodefaults__')
            nodefaults = vals{end};
            vals(end-1:end) = [];
        end
        
        % check for skip flag (skip leading skippable args in positional argument lists)
        if length(vals)>1 && isequal(vals{end-1},'__arg_skip__')
            skip = vals{end};
            vals(end-1:end) = [];
        end

    else
        report_type = 'none';
    end
        
    % if Format is a function, run it to reformat vals and fmt
    if isa(fmt,'function_handle')
        switch nargin(fmt)
            case 1
                if nargout(fmt) == 1
                    vals = feval(fmt,vals); 
                    fmt = 0;
                else
                    [vals,fmt] = feval(fmt,vals);
                end
            case 2
                % first expand the spec
                tmpspec = hlp_microcache('spec',@expand_spec,compressed_spec,'lean',~nodefaults,hlp_getcaller(2));
                % then call the function with it
                if nargout(fmt) == 1
                    vals = feval(fmt,vals,tmpspec);
                    fmt = 0;
                else
                    [vals,fmt] = feval(fmt,vals,tmpspec);
                end
            otherwise
                error('The given formatting function expects an unsupported number of inputs (only 1 or 2 inputs supported).');
        end            
        if isa(fmt,'function_handle')
            error('If a given formatting function returns a new format description as second output, that description must be numeric (e.g., [0 Inf] or 0).'); end
    end
    
    % get places where structs occur in the Values
    structmask = cellfun('isclass',vals,'struct');
end


% check if arg_define is being called in direct mode
function direct_mode = is_direct_mode(vals,structmask)
    if length(vals)>1 && isequal(vals(end-1:end),{'arg_direct',true})
        direct_mode = true;
    elseif ~isempty(vals) && isfield(vals{end},'arg_direct')
        direct_mode = vals{end}.arg_direct;
    else
        direct_mode = false;
        % if a report is requested, we cannot be in direct mode
        indices = find(structmask | strcmp(vals,'arg_direct'));
        for k = indices(end:-1:1)
            if ischar(vals{k}) && k<length(vals)
                % found it in the NVPs
                direct_mode = vals{k+1};
                break;
            elseif isfield(vals{k},'arg_direct')
                % found it in a struct
                direct_mode = vals{k}.arg_direct;
                break;
            end
        end
    end
end


% directly assign names to values
function [nvps,outstruct] = assign_direct(vals,structmask,make_struct)
    if make_struct
        % the variables are returned in the struct named result
        if ~isscalar(structmask)
            % obtain flat NVP list
            nvps = flatten_structs(vals,structmask);
            % split names/values and find indices of the last assignment for each name
            nvps = reshape(nvps,2,[]);
            [s,indices] = sort(nvps(1,:));
            indices(strcmp(s(1:end-1),s(2:end))) = [];
            % build & return the struct
            outstruct = cell2struct(nvps(2,indices),nvps(1,indices),2);
        else
            % a single struct was passed in: pass it right through
            nvps = vals;
            outstruct = nvps{1};
        end
    else
        % no output argument: generate NVP list for subsequent assignment to function workspace
        nvps = flatten_structs(vals,structmask);
        outstruct = [];
    end
end


% cached wrapper around process_spec
function varargout = process_spec_cached(caller_name,spec,report_type,assign_defaults,perform_namecheck)
    persistent cache;
    key = [report_type 'a'+assign_defaults];
    try
        % try to load the cached result for the given caller
        result = cache.(caller_name).(key);
        % return if it matches the current input
        if isequal(result.spec,spec) || isequalwithequalnans(result.spec,spec) %#ok<FPARK>
            varargout = result.output; 
            return; 
        end
    catch %#ok<CTCH>
    end
    % otherwise fall back to microcache and overwrite
    [varargout{1:nargout}] = hlp_microcache('spec',@process_spec,spec,report_type,assign_defaults,perform_namecheck);
    cache.(caller_name).(key) = struct('output',{varargout},'spec',{spec});
end


% expand a specification from a cell array into a struct array
function [spec,flat_names,first_names,name2idx,leading_skippable,checks] = process_spec(compressed_spec,report_type,assign_defaults,perform_namecheck)
    caller_name = hlp_getcaller(2);
    
    % expand the compressed cell-array spec into a full-blown struct-array spec
    spec = expand_spec(compressed_spec,'rich',assign_defaults,caller_name);
    
    % obtain the argument names and a flat list thereof
    arg_names = {spec.names};
    flat_names = [arg_names{:}];
    first_names = cellfun(@(n)n{1},arg_names,'UniformOutput',false);

    % create a name/index remapping table
    name2idx = struct();
    for n=1:length(arg_names)
        for k=1:length(arg_names{n})
            name2idx.(arg_names{n}{k}) = n; end
    end
    
    % find the number of leading skippable arguments
    leading_skippable = find(~[spec.skippable false],1)-1;

    % check for duplicate argument names in the Specification
    sorted_names = sort(flat_names);
    duplicates = unique(sorted_names(strcmp(sorted_names(1:end-1),sorted_names(2:end))));
    if ~isempty(duplicates)
        error(['The names ' hlp_tostring(duplicates) ' refer to multiple arguments.']); end
    
    % if required, check for name clashes with functions on the path
    % (this is due to a deficiency in MATLAB's handling of variables that were assigned to a 
    % function's scope from a subfunction, which are prone to clashes with functions on the path...)
    if perform_namecheck && strcmp(report_type,'none')
        try
            validate_varnames(caller_name,first_names);
        catch e
            disp_once('The function validate_names failed; reason: %s',e.message);
        end
    end
    
    % determine which checks are necessary
    checks.mandatory = check_property(spec,'value',mandatory);
    checks.unassigned = check_property(spec,'value',unassigned);
    checks.expression = check_property(spec,'type','expression');
    checks.conversion = check_property(spec,'to_double',true);
end


% expand the given spec cell array into a struct array
function spec = expand_spec(spec,report_type,assign_defaults,caller_name)
    % evaluate the cells into specifier structs
    if all(cellfun('isclass',spec,'cell'))
        spec = cellfun(@(s)feval(s{1},report_type,s{2}{:}),spec); end
    % make sure that spec has the correct fields, even if empty
    if isempty(spec)
        spec = arg_specifier;
        spec = spec([]);
    end
    % optionally assign the default values
    if assign_defaults
        for s=1:length(spec)
            for def=spec(s).defaults
                spec(s) = assign_value(spec(s),def{1},report_type,caller_name,false,false); end
        end
    end    
end


% recursively check whether the given property is set to the given value in any of the spec entries
function ismatch = check_property(spec,name,value)
    if isempty(spec)
        ismatch = false;
    else
        if ischar(value)
            ismatch = any(strcmp({spec.(name)},value));
        elseif isscalar(value) && islogical(value)
            ismatch = any([spec.(name)]==value);
        else
            ismatch = any(cellfun(@(x)isequal(x,value),{spec.(name)}));
        end
        if ~ismatch
            for k=find(~cellfun('isempty',{spec.children}))
                ismatch = check_property(spec(k).children,name,value);
                if ismatch
                    return; end
            end
            for k=find(~cellfun('isempty',{spec.alternatives}))
                ismatch = check_property([spec(k).alternatives{:}],name,value);
                if ismatch
                    return; end            
            end
        end
    end
end


% transform the cell array of input arguments (vals) to a pure list of name-value pairs (NVPs)
function nvps = arguments_to_nvps(caller,fmt,vals,structmask,flat_names,first_names,skipped_positionals)
    permitted_names = [flat_names {'arg_selection','arg_direct'}];
    
    % get the call signature; this is vals with everything replaced by [] that is neither a valid argument 
    % name nor a struct, and structs replaced by a sorted cell array of their field names
    charmask = cellfun('isclass',vals,'char');
    signature = vals;
    signature(~(charmask|structmask)) = {[]};
    for k=find(charmask)
        if ~any(strcmp(signature{k},permitted_names))
            signature{k} = sprintf('%u',k); end
    end
    if any(structmask)
        signature(structmask) = cellfun(@fieldnames,signature(structmask),'UniformOutput',false); end
    
    % determine the number n of arguments that were specified positionally
    [n,violations,ignored] = hlp_nanocache(caller,10,@get_num_positionals,fmt,length(first_names),signature,structmask,permitted_names,skipped_positionals);
    
    % emit errors and/or diagnostic warnings
    if isnan(n)
        for k=find(cellfun(@(s)all(s>='0'&s<='9'),violations))
            violations{k} = vals{str2num(violations{k})}; end %#ok<ST2NM>
        error([hlp_getcaller(2) ':arg_define:invalid_arguments'],['Some of the specified arguments do not appear in the argument specification; ' hlp_tostring(violations) '.']);
    elseif ~isempty(ignored)
        for k=find(cellfun(@(s)all(s>='0'&s<='9'),violations))
            violations{k} = vals{str2num(violations{k})}; end %#ok<ST2NM>
        caller_name = hlp_getcaller(2);
        warn_once([caller_name ':arg_define:possible_conflict'],'arg_define() in %s: Possible parameter conflict -- both unrecognized parameters %s and matching names %s passed in. Assuming that the function is called with %u positional arguments. This warning will not be repeated for this MATLAB session.',caller_name,hlp_tostring(violations),hlp_tostring(ignored),n);
    end
    
    % generate a flat list of name-value pairs
    if n
        % optionally pad vals with implicit []'s for each positional to skip
        if skipped_positionals
            vals = [repmat({[]},1,skipped_positionals) vals]; end
        % prepend the positionals and their associated names
        nvps = [reshape([first_names(1:n);vals(1:n)],1,[]) flatten_structs(vals(n+1:end),structmask(n+1:end))];
    elseif any(structmask)
        nvps = flatten_structs(vals,structmask);
    else
        nvps = vals;
    end
end


% determine the number of positional arguments in a given call signature
function [num_positionals,violations,ignored] = get_num_positionals(fmt,spec_length,signature,structmask,permitted_names,skipped_positionals)
    % preprocess fmt
    max_positionals = min(spec_length,length(signature)+skipped_positionals);
    if isempty(fmt)
        % empty Format means that any number of positionals are permitted
        fmt = 0:max_positionals;
    else
        % make sure that 0 positional arguments are always permitted
        if ~any(fmt==0)
            fmt = [0 fmt]; end
        % remap Inf to the max # of positionals
        fmt(fmt==Inf) = max_positionals;
        % remove anything else that is out of bounds
        fmt(fmt>max_positionals | fmt<0) = [];
        % make sure that it is sorted and contains no duplicates
        fmt = unique(sort(fmt));
    end

    % replace the cell arrays in the signature by dummy structs, so we can use flatten_structs
    for k=find(structmask)
        signature{k} = cell2struct(cell(1,length(signature{k})),signature{k},2); end
    
    % for each permitted number of positional arguments n in fmt...
    num_positionals = NaN;
    violations = {};
    ignored = {};
	for n = fmt %#ok<*ALIGN>
        if skipped_positionals && length(fmt)>1 && n==fmt(2)
            % pad the signature with dummy positionals if num_skippable is nonzero
            signature = [repmat({[]},1,skipped_positionals) signature]; %#ok<AGROW>
            structmask = [false(1,skipped_positionals) structmask]; %#ok<AGROW>
        end
        % interpret the arguments after the first n positions als NVPs/structs
		nvps = flatten_structs(signature(n+1:end),structmask(n+1:end));
        % perform type check
		if iscellstr(nvps(1:2:end))
		    % make sure that only admissible names are used
		    inadmissible = fast_setdiff(nvps(1:2:end),permitted_names);
   		    if isempty(inadmissible)
                % NVP arguments are all admissible: use this as num_positionals
                num_positionals = n;                
                if n>0
                    % assuming n positional arguments; collect any admissible names among those, i.e.,
                    % names that would have been valid argument names but were ignored
                    fmt(fmt>=n) = [];
                    for k=fmt
                        nvps = flatten_structs(signature(k+1:n),structmask(k+1:n));
                        if iscellstr(nvps(1:2:end))
                            ignored = [ignored intersect(nvps(1:2:end),permitted_names)]; %#ok<AGROW>
                        else
                            ignored = [ignored intersect(nvps(cellfun('isclass',nvps,'char')),permitted_names)]; %#ok<AGROW>
                        end
                    end
                    ignored = unique(ignored);
                end                
                violations = unique(violations);
                return;
            else
                % inadmissible names: collect diagnostic information
                violations = [violations inadmissible]; %#ok<AGROW>
		    end
        else
            violations = [violations {'invalid sequence of names and structs'}]; %#ok<AGROW>
        end
    end
    violations = unique(violations);
end


% assign the values in an NVP list to the spec
function spec = assign_nvps(spec,nvps,name2idx,report_type,caller_name,deprecation_warning)
    if ~strcmp(report_type,'rich')
        report_type = 'lean'; end
    % note: this part needs to be changed to accommodate arg_ref()
    for k=1:2:length(nvps)
        try 
            idx = name2idx.(nvps{k});
            spec(idx) = assign_value(spec(idx),nvps{k+1},report_type,caller_name,true,deprecation_warning);
        catch e
            if ~strcmp(e.identifier,'MATLAB:nonExistentField')
                rethrow(e); end
            % this is an internal sanity check; if this message is triggered some internal error has
            % occurred
            if any(strcmp({spec.first_name},nvps{k}))
                error(['Cannot insert a duplicate field into the specification: ' nvps{k}]); end
            % append it to the spec (note: this could use some optimization... it would be better
            % if the spec automatically contained the arg_selection field)
            spec(end+1) = cached_argument(nvps{k},nvps{k+1});
            name2idx.(nvps{k}) = length(spec);
        end
    end
end


% assign a given value to a single element in a spec
function spec = assign_value(spec,newvalue,report_type,caller_name,nodefaults,deprecation_warning)
    skip_arg = {'__arg_skip__',true};
    nodefault_arg = {'__arg_nodefaults__',true};
    % check whether this value is assignable
    if ~isequal(newvalue,'__arg_unassigned__') && ~(~spec.empty_overwrites && (isempty(newvalue) || isequal(newvalue,'__arg_mandatory__')))
        % warn about deprecation
        if deprecation_warning && spec.deprecated && ~isequal_weak(spec.value,newvalue)
            warn_deprecation(spec,caller_name); end
        % perform assignment
        if isempty(spec.mapper)
            spec.value = newvalue;
        else
            % apply the mapper to get the selection key and the value pack
            [key,value] = spec.mapper(newvalue,spec.range,spec.names);
            spec.value = key;
            if isempty(key)
                % arg_sub
                pos = 1;
                source_fields = spec.children;
            elseif islogical(key) || isnumeric(key)
                % arg_subtoggle
                pos = key+1;
                source_fields = spec.alternatives{pos};
            elseif ischar(key)
                % arg_subswitch
                pos = find(strcmp(key,spec.range));
                source_fields = spec.alternatives{pos};
            else
                error('Unsupported mapper key class.');
            end                        

            % parse the value into a spec struct array
            if isequal(key,false)
                spec.value = false;
                spec.children = cached_argument('arg_selection',false);
            else
                if nodefaults
                    % parse just the value, without applying defaults
                    value = arg_report(report_type,spec.sources{pos},[value skip_arg nodefault_arg]);
                    % selectively override the current source fields with the value
                    spec.children = override_fields(source_fields,value);
                else
                    % replace the children by the result
                    spec.children = arg_report(report_type,spec.sources{pos},[value skip_arg]);
                end
                % override flags
                spec.children = override_flags(spec.children,spec.reflag{pos}{:});
                if ~isempty(key)
                    % set/append selector child
                    selection_arg = strcmp('arg_selection',{spec.children.first_name});
                    if any(selection_arg)
                        spec.children(selection_arg).value = key;
                    else
                        spec.children = [spec.children,cached_argument('arg_selection',key)];
                    end
                    % update alternatives
                    spec.alternatives{pos} = spec.children; 
                end
            end
        end
    end
end


% emit a deprecation warning
function warn_deprecation(spec,caller_name)
    if iscell(spec.help)
        if length(spec.help) == 1
            help = spec.help{1};
        else
            help = [spec.help{1} '.' spec.help{2}];
        end
    else
        help = spec.help;
    end
    if length(spec.names) > 1
        name = [spec.names{1} '/' spec.names{2}];
    else
        name = spec.names{1};
    end
    disp_once(['Using deprecated argument "' name '" in function ' caller_name ' (help: ' help ').']); 
end