function spec = expand_argsubtoggle(varargin)
% Internal: expand the output of an arg_subtoggle(...) declaration function into an argument specifier.
% Specifier = expand_argsubtoggle(ReportType, ...)
%
% The argument declaration functions used in an arg_define declaration produce a compact and
% low-overhead representation that first needs to be expanded into a full specifier struct that
% contains all properties of the argument and can be processed by other functions. This function
% performs the expansion for the arg_subtoggle() arguments.
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

    % the result of build_specifier is cached in memory for efficient processing of repeated calls
    spec = hlp_microcache('arg',@build_specifier,varargin{:});
end


% expand the arg_subtoggle(...) declaration line into an argument specifier (see expand_argsub for a simpler example)
function spec = build_specifier(reptype,names,defaults,source,help,varargin)
    % set defaults
    if nargin < 4 || isempty(source)
        error('BCILAB:args:no_options','The Source argument for arg_subtoggle() may not be omitted.'); end 
    if nargin < 5
        help = ''; end
    [fmt,reflag,suppress,permit_positionals,skip_noreps,alternative_defaults] = deal([0 Inf],{},{},true,true,{});
    
    % extract special options
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'fmt','reflag','suppress','permit_positionals','skip_noreps','alternative_defaults'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end
    
    % handle the 'suppress' option (by appending to reflag)
    for n=suppress
        reflag = [reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % initialize the specification struct
    spec = arg_specifier('head',@arg_subtoggle,'names',names,'help',help,'mapper',@map_argsubtoggle,varargin{:},'type','logical','shape','scalar');
    
    % parse the 'mapper'
    if nargin(spec.mapper) == 1
        spec.mapper = @(x,y,z) spec.mapper(x); end
    
    % parse the Source
    source = generate_source(fmt,source);
    
    % rewrite 'skip_noreps' and Defaults into more convenient forms
    skip_noreps = quickif(skip_noreps,{'__arg_skip__',true},{});
    [default_sel,default_val] = spec.mapper(defaults);
    if ~isempty(default_val)
        default_spec = arg_report('rich',source,default_val);
        default_val = arg_tovals(default_spec,[],'cell');
    else
        default_spec = [];
    end
    
    % set up the assigner
    spec.assigner = @(spec,value) assign_argsubtoggle(spec,value,reptype,source,default_spec,default_val,reflag,permit_positionals,skip_noreps);
    
    % populate the alternatives in case of a rich spec (rich spec contains alternative branches)
    if strcmp(reptype,'rich')
        if default_sel
            spec.alternatives{1} = override_flags(cached_selector(false),reflag{:});
        else
            spec.alternatives{2} = override_flags([arg_report(reptype,source,alternative_defaults) cached_selector(true)],reflag{:});
        end
    end
    
    % assign the default
    spec = assign_argsubtoggle(spec,defaults,reptype,source,[],{},reflag,permit_positionals,skip_noreps);
end


% this function maps an argument list onto a binary flag (enabled status) and the cell array of 
% name-value pairs / structs to assign
function [selected,args] = map_argsubtoggle(args)
    if isequal(args,'on') || isequal(args,{})
        selected = true;
        args = {};
    elseif isequal(args,'off') || (isempty(args) && isnumeric(args))
        selected = false;
        args = {};
    elseif isfield(args,'arg_selection') && isscalar(args)
        selected = args.arg_selection;
        args = {args};
    elseif iscell(args)        
        if isscalar(args) && isfield(args{1},'arg_selection') && isscalar(args{1}.arg_selection)
            selected = args{1}.arg_selection;
        elseif isequal(args,{'arg_selection',false})
            selected = false;
            args = {};
        elseif isequal(args,{'arg_selection',true})
            selected = true;
            args = {};
        else
            pos = find(strcmp('arg_selection',args(1:end-1)),1,'last');
            if isempty(pos)
                selected = true;
            else
                [selected,args] = deal(args{pos+1},args([1:pos-1 pos+2:end]));
            end
        end
    else
        selected = true;
        args = {args};
    end
end


% function used to assign a value to the argument
function spec = assign_argsubtoggle(spec,value,reptype,source,default_spec,default_val,reflag,permit_positionals,skip_noreps)
    % skip unassignable values
    if isequal(value,'__arg_unassigned__') || (~spec.empty_overwrites && (isempty(value) || isequal(value,'__arg_mandatory__')))
        return; end
    
    % apply the mapper
    [spec.value,value] = spec.mapper(value);
    
    % assign the children
    if spec.value
        if ~isempty(default_val) && permit_positionals
            % parse the values into a struct and retain only the difference from the rich default_val
            diffvalue = arg_tovals(arg_diff(default_spec,arg_report('lean',source,[value skip_noreps])),[],'cell');
            % now parse the default_val with values partially overriding and assign result to children (optionally reflagged)
            spec.children = arg_report(reptype,source,[default_val,diffvalue]);
        else
            % optimization: can just concatenate default_val and value
            spec.children = arg_report(reptype,source,[default_val value skip_noreps]);
        end
        % set or append selector
        selection_arg = strcmp('arg_selection',{spec.children.first_name});
        if any(selection_arg)
            spec.children(selection_arg).value = spec.value;
        else
            spec.children = [spec.children,cached_selector(spec.value)];
        end
    else
        spec.children = cached_selector(spec.value);
    end
    
    % override flags
    spec.children = override_flags(spec.children,reflag{:});
    
    % also override the corresponding entry in alternatives
    spec.alternatives{1+spec.value} = spec.children;    
end


% returns a cached selector argument specifier
function result = cached_selector(selected)
    persistent cache;
    try
        result = cache.(char('a'+selected));
    catch %#ok<CTCH>
        result = arg_nogui('arg_selection',selected);
        result = feval(result{1},[],result{2}{:});        
        cache.(char('a'+selected)) = result;
    end
end
