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
    if ischar(suppress)
        suppress = {suppress}; end
    for n=suppress
        reflag = [reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % initialize the specification struct
    spec = arg_specifier('head',@arg_subtoggle,'names',names,'help',help,'mapper',@map_argsubtoggle,varargin{:},'type','logical','shape','scalar');
    
    % parse the 'mapper'
    if nargin(spec.mapper) == 1
        spec.mapper = @(x,y,z) spec.mapper(x); end
    
    % parse the Source
    source = generate_source(fmt,source);
    
    % rewrite 'skip_noreps'
    skip_noreps = quickif(skip_noreps,{'__arg_skip__',true},{});
    
    % set up the assigner
    spec.assigner = @(spec,value) assign_argsubtoggle(spec,value,reptype,source,reflag,permit_positionals,skip_noreps);
    
    % initialize the contents
    spec.contents = repmat({{}},1,2);
    
    % populate the alternatives in case of a rich spec (rich spec contains alternative branches)
    if strcmp(reptype,'rich')
        spec.alternatives = cell(1,2);
        spec = spec.assigner(spec,[]);
        spec = spec.assigner(spec,alternative_defaults);
    end
    
    % assign the defaults
    spec = spec.assigner(spec,defaults);
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
function spec = assign_argsubtoggle(spec,value,reptype,source,reflag,permit_positionals,skip_noreps)
    % skip unassignable values
    if isequal(value,'__arg_unassigned__') || (~spec.empty_overwrites && (isempty(value) || isequal(value,'__arg_mandatory__')))
        return; end
    
    % apply the mapper
    [spec.value,value] = spec.mapper(value);
    idx = spec.value+1;
    
    if spec.value
        value = hlp_microcache('argparse',@arg_report,'parse',source,[value skip_noreps]);
        if length(spec.alternatives)>=idx && length(spec.alternatives{idx})>1 
            diffvalue = arg_diff(spec.alternatives{idx},value);
        else
            diffvalue = value;
        end
        diffvalue = arg_tovals(diffvalue,[],'cell',false,false,false,false);
        spec.contents{idx} = [spec.contents{idx} diffvalue];
        spec.children = arg_report(reptype,source,spec.contents{idx});    

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
    
    % override the corresponding entry in alternatives
    spec.alternatives{idx} = spec.children;
end


% returns a cached selector argument specifier
function result = cached_selector(selected)
    persistent cache;
    try
        result = cache.(char('a'+selected));
    catch %#ok<CTCH>
        result = arg_nogui('arg_selection',selected,[],[],'assigned',true);
        result = feval(result{1},[],result{2}{:});
        cache.(char('a'+selected)) = result;
    end
end
