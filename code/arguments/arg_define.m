function outstruct = arg_define(vals,varargin)
% Declare function arguments with support for introspection.
% Struct = arg_define(Values, Specification...)
% Struct = arg_define(Format, Values, Specification...)
%
% This is a replacement for the parameter declaration line of a function (named the "Function" in
% the following). The content of Values (a cell array of values, typically the list of inputs
% received by the Function), is assigned to the fields in the output Struct, or to variables in the
% Function's workspace. Parsing is done according to a Specification of argument names and their
% order (optionally obeying a custom Format description).
%
% Values can be a list of a fixed number of positional arguments (i.e., the typical MATLAB calling
% format), optionally followed by a list of name-value pairs (NVPs, e.g., as the format accepted by
% figure()). The name-value pairs may be interleaved with structs, i.e., one may pass a mix of
% 'name',value,STRUCT,'name',value,'name',value, ...). Note that this is not ambiguous. Passing in
% everything as NVPs/structs is generally allowed. Only names that are listed in the Specification
% may be used in the NVPs.
%
% The difference to other NVP parsing functions is that arguments defined via arg_define() can be
% reported to outside functions (queried using arg_report()). The resulting argument specification
% struct can be used to draw a GUI, store settings, generate help text, or be processed otherwise.
%
% In:
%   Format : Note that this argument can be omitted. Optionally this is the number of allowed
%            positional arguments of the Function (default: [0 Inf]), i.e., the number of leading
%            arguments that can be passed to the Function as positional arguments, while the
%            remaining arguments are interpreted as a list of name-value pairs and/or structs.
%
%            Special uses: 
%            * If Format is a vector of multiple values, then the smallest number for which the 
%              Function's inputs are consistent, in the sense that each name in Values matches a
%              name in the Specification, is accepted.
%            * If Format is numeric but does not contain a 0, a 0 will be implicitly prepended.
%            * If Format is [], any number of positional arguments is allowed (i.e., Format defaults
%              to 0:length(Specification)), although this is discouraged in practice as it may
%              result in mis-parsing when mis-spelled arguments are passed in (arg_define will then
%              assume that anything up to the mis-spelled argument is passed in positionally).
%            * If Format is a function handle, the given function can be used to transform the
%              Values prior to any other processing into a new Values cell array. The function may
%              optionally specify a new (numeric) Format as its second output argument (if not
%              specified, this is 0). The function may optionally take in the argument specification
%              (struct array of arg_specifier), although some details of this data structure are
%              internal and subject to change in future releases (particularly fields marked as
%              INTERNAL).
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
% Technical Note:
%   The Function may be called with the request to deliver the parameter specification as opposed to
%   following the normal execution flow (this is done by arg_report). arg_define responds to this
%   task by throwing an exception of the type 'BCILAB:arg:report_args' using arg_issuereport, which
%   is caught by arg_report to obtain the specification struct.
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

    % first parse the inputs
    [fmt,vals,structmask,spec,report,skip] = process_inputs(vals,varargin);
    
    % check if we are in direct mode (fast shortcut)
    if strcmp(report,'none') && is_direct_mode(vals,structmask)
        % in direct mode we assign only those variables that were passed in
        [nvps,outstruct] = assign_direct(vals,structmask,nargout);
    else
        % otherwise we perform full parsing
        
        % get the name of the calling function (for diagnostics and cache lookups)
        caller = hlp_getcaller(2);
        if ~isvarname(caller)
            caller = 'anonymous'; end
        
        % process the specification into a struct array and extract some properties
        [spec,flat_names,first_names,name2idx,leading_skippable,checks] = process_spec_cached(caller,spec,report,nargout==0);

        % convert vals to a canonical list of name-value pairs (NVPs)
        nvps = arguments_to_nvps(caller,fmt,vals,structmask,flat_names,first_names,skip*leading_skippable);

        % assign the NVPs to the spec
        spec = assign_values(spec,nvps,name2idx);

        % generate outputs from spec
        switch report
            case 'none'
                outstruct = arg_tovals(spec,[],'struct',checks.mandatory,checks.unassigned,checks.expression,checks.conversion);
                nvps = reshape([fieldnames(outstruct)';struct2cell(outstruct)'],1,[]);
            case 'vals'
                arg_issuereport(arg_tovals(spec,[],'struct',checks.mandatory,checks.unassigned,checks.expression,checks.conversion));
            case 'nvps'
                arg_issuereport(arg_tovals(spec,[],'cell',checks.mandatory,checks.unassigned,checks.expression,checks.conversion));
            case {'lean','rich'}
                arg_issuereport(spec);
            otherwise
                error(['Unrecognized report type requested: ' report]);
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
function [fmt,vals,structmask,spec,report,skip] = process_inputs(vals,spec)
    if iscell(vals)
        % no Format specifier was given: use default
        fmt = [0 Inf];
    else
        % a Format specifier was given as first argument (need to shift remaining arguments by 1)
        fmt = vals;
        vals = spec{1};
        spec(1) = [];
    end
    
    % extract report type
    if length(vals)>1 && isequal(vals{end-1},'__arg_report__')
        report = vals{end};
        % check for some report types that can be handled immediately
        switch report
            case 'raw'
                arg_issuereport(spec);
            case 'properties'
                arg_issuereport(struct());
            case {'lean','rich'}
                if length(vals) == 2
                    arg_issuereport(hlp_microcache('spec',@expand_spec,spec,report)); end
            case 'handle'
                error('To make function handles accessible, use the function expose_handles().');
        end
        vals(end-1:end) = [];
    else
        report = 'none';
    end
    
    % extract skip flag (skip skippable positional arguments)
    if length(vals)>1 && isequal(vals{end-1},'__arg_skip__')
        skip = vals{end};
        vals(end-1:end) = [];
    else
        skip = false;
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
                tmpspec = hlp_microcache('spec',@expand_spec,spec,'lean');
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
function varargout = process_spec_cached(caller,spec,type,perform_namecheck)
    persistent cache;
    try
        % try to load the cached result for the given caller
        result = cache.(caller).(type);
        % return if it matches the current input
        if isequal(result.spec,spec) || isequalwithequalnans(result.spec,spec)
            varargout = result.output; 
            return; 
        end
    catch %#ok<CTCH>
    end
    % otherwise fall back to microcache and overwrite
    [varargout{1:nargout}] = hlp_microcache('spec',@process_spec,spec,type,perform_namecheck);
    cache.(caller).(type) = struct('output',{varargout},'spec',{spec});
end


% expand a specification from a cell array into a struct array
function [spec,flat_names,first_names,name2idx,leading_skippable,checks] = process_spec(spec,report,perform_namecheck)
    % first expand the cell-array spec into a struct-array spec
    spec = expand_spec(spec,quickif(strcmp(report,'rich'),'rich','lean'));

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
    % (this is due to a deficiency in MATLAB's handling of variables that were assigned to a function's scope
    % from the outside, which are prone to clashes with functions on the path...)
    if perform_namecheck && strcmp(report,'none')
        try
            validate_names(first_names);
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


% check for name clashes (once)
function validate_names(varnames)
    persistent already_checked;
    if isempty(already_checked)
        already_checked = {}; end
    for n = fast_setdiff(varnames,already_checked)
        curname = n{1};
        already_checked{end+1} = curname; %#ok<AGROW>
        existing_func = which(curname);
        if ~isempty(existing_func)
            if isempty(strfind(existing_func,'Java method'))
                [path_part,file_part,ext_part] = fileparts(existing_func);
                if ~any(strncmp('@',hlp_split(path_part,filesep),1))
                    % If this happens, it means that there is a function in one of the directories in
                    % MATLAB's path which has the same name as an argument of the specification. If this
                    % argument variable is copied into the function's workspace by arg_define, most MATLAB
                    % versions will (incorrectly) try to call that function instead of accessing the
                    % variable. I hope that they handle this issue at some point. One workaround is to use
                    % a longer argument name (that is less likely to clash) and, if it should still be
                    % usable for parameter passing, to retain the old name as a secondary or ternary
                    % argument name (using a cell array of names in arg()). The only really good
                    % solution at this point is to generally assign the output of arg_define to a
                    % struct.
                    disp([hlp_getcaller(4) ': The argument name "' curname '" clashes with the function "' [file_part ext_part] '" in directory "' path_part '"; it is strongly recommended that you either rename the function or remove it from the path.']);
                end
            else
                % these Java methods are probably spurious "false positives" of the which() function
                disp([hlp_getcaller(4) ': There is a Java method named "' curname '" on your path; if you experience any name clash with it, please report this issue.']);
            end
        end
    end
end


% substitute any structs in place of a name-value pair into the name-value list
function args = flatten_structs(args,structmask)
    if any(structmask)
        persistent cache; %#ok<TLEV>
        try
            % try to look up splice points from cache
            if length(structmask) < 62
                field = char('a'+structmask);
            else
                str = java.lang.String(char('a'+structmask));
                field = ['a' num2str(str.hashCode()+2^31)];
            end
            splicepos = cache.(field);
        catch %#ok<CTCH>
            % pre-calculate splice points from structmask
            % and cache results
            splicepos = [];
            k = 1;
            while k <= length(args)
                if isstruct(args{k})
                    splicepos(end+1) = k; %#ok<AGROW>
                    k = k+1; % struct case
                else
                    k = k+2; % NVP case
                end
            end
            splicepos = splicepos(end:-1:1);
            cache.(field) = splicepos;
        end

        % splice structs in
        for k = splicepos
            args = [args(1:k-1) reshape([fieldnames(args{k}) struct2cell(args{k})]',1,[]) args(k+1:end)]; end
    end
end


% expand the given spec cell array into a struct array
function spec = expand_spec(spec,report)
    % evaluate the cells into specifier structs
    if all(cellfun('isclass',spec,'cell'))
        spec = cellfun(@(s)feval(s{1},report,s{2}{:}),spec); end
    % make sure that spec has the correct fields, even if empty
    if isempty(spec)
        spec = arg_specifier;
        spec = spec([]);
    end
end


% recursively check whether the given property matches the given value in any of the spec entries
function ismatch = check_property(spec,name,value)
    if isempty(spec)
        ismatch = false;
    else
        ismatch = any(cellfun(@(x)isequal(x,value),{spec.(name)}));
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
        signature(structmask) = cellfun(@(s)fieldnames(s),signature(structmask),'UniformOutput',false); end
    
    % determine the number n of arguments that were specified positionally
    [n,violations,ignored] = hlp_nanocache(caller,10,@get_num_positionals,fmt,length(first_names),signature,structmask,permitted_names,skipped_positionals);
    
    % emit errors and/or diagnostic warnings
    if isnan(n)
        for k=find(cellfun(@(s)all(s>='0'&s<='9'),violations))
            violations{k} = vals{str2num(violations{k})}; end %#ok<ST2NM>
        error([hlp_getcaller(3) ':arg_define:invalid_arguments'],['Some of the specified arguments do not appear in the argument specification; ' hlp_tostring(violations) '.']);
    elseif ~isempty(ignored)
        for k=find(cellfun(@(s)all(s>='0'&s<='9'),violations))
            violations{k} = vals{str2num(violations{k})}; end %#ok<ST2NM>
        caller_name = hlp_getcaller(3);
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
function spec = assign_values(spec,nvps,name2idx)
    % note: this part needs to be changed to accommodate arg_ref()
    for k=1:2:length(nvps)
        try
            idx = name2idx.(nvps{k});
            newvalue = nvps{k+1};
            if spec(idx).deprecated && ~isequal_weak(spec(idx).value,newvalue)
                disp_once(['Using deprecated argument "' nvps{k} '" in function ' hlp_getcaller(2) ' (help: ' hlp_tostring(spec(idx).help) ').']); end
            if ~isempty(spec(idx).assigner) 
                spec(idx) = spec(idx).assigner(spec(idx),newvalue);
            elseif ~isequal(newvalue,'__arg_unassigned__') && ~(~spec(idx).empty_overwrites && (isempty(newvalue) || isequal(newvalue,'__arg_mandatory__')))
                spec(idx).value = newvalue;
            end
        catch e
            if ~strcmp(e.identifier,'MATLAB:nonExistentField')
                rethrow(e); end
            % this is an internal sanity check; if this message is triggered some internal error has
            % occurred
            if any(strcmp({spec.first_name},nvps{k}))
                error(['Cannot insert a duplicate field into the specification: ' nvps{k}]); end
            % append it to the spec (note: this might need some optimization... it would be better
            % if the spec automatically contained the arg_selection field)
            tmp = arg_nogui(nvps{k},nvps{k+1});
            spec(end+1) = feval(tmp{1},[],tmp{2}{:}); %#ok<AGROW>
        end
    end
end
