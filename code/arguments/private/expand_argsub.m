function spec = expand_argsub(varargin)
% Internal: expand the output of an arg_sub(...) declaration function into an argument specifier.
% Specifier = expand_argsub(ReportType, ...)
%
% The argument declaration functions used in an arg_define declaration produce a compact and
% low-overhead representation that first needs to be expanded into a full specifier struct that
% contains all properties of the argument and can be processed by other functions. This function
% performs the expansion for the arg_sub() arguments.
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


% expand the arg_sub(...) declaration line into an argument specifier (see expand_arg for a simpler example)
function spec = build_specifier(reptype,names,defaults,source,help,varargin)
    % set defaults
    if nargin < 2
        names = {}; end
    if nargin < 3
        defaults = {}; end
    if nargin < 4
        source = {}; end
    if nargin < 5
        help = ''; end    
    [fmt,reflag,suppress,permit_positionals,skip_noreps] = deal([0 Inf],{},{},true,true);
    
    % extract special options
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'fmt','reflag','suppress','permit_positionals','skip_noreps'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end    
    
    % handle the 'suppress' option (by appending to reflag)
    for n=suppress
        reflag = [reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % initialize the specification struct
    spec = arg_specifier('head',@arg_sub,'names',names,'help',help,varargin{:},'value',[],'type','char','shape','row');
    
    % parse the Source
    source = generate_source(fmt,source);
    
    % rewrite 'skip_noreps' and Defaults into more convenient forms
    skip_noreps = quickif(skip_noreps,{'__arg_skip__',true},{});
    if ~isempty(defaults)
        default_spec = arg_report('rich',source,defaults);
        default_val = arg_tovals(default_spec,[],'cell');
    else
        default_spec = [];
        default_val = {};
    end
    
    % set up the assigner
    spec.assigner = @(spec,value) assign_argsub(spec,value,reptype,source,default_spec,default_val,reflag,permit_positionals,skip_noreps);
    
    % assign the default
    spec = assign_argsub(spec,default_val,reptype,source,[],{},reflag,permit_positionals,skip_noreps);
end


% function to perform the value assignment
function spec = assign_argsub(spec,value,reptype,source,default_spec,default_val,reflag,permit_positionals,skip_noreps)
    % skip unassignable values
    if isequal(value,'__arg_unassigned__') || (~spec.empty_overwrites && (isempty(value) || isequal(value,'__arg_mandatory__')))
        return; end
    % make sure that value is a cell array (otherwise we cannot concatenate it below)
    if ~iscell(value)
        value = {value}; end
    
    if ~isempty(default_val) && permit_positionals
        % parse the values into a struct and retain only the difference from the rich default_val
        diffvalue = arg_tovals(arg_diff(default_spec,arg_report('lean',source,[value skip_noreps])),[],'cell');
        % now parse the default_val with values partially overriding and assign result to children (optionally reflagged)
        spec.children = override_flags(arg_report(reptype,source,[default_val diffvalue]),reflag{:});
    else
        % optimization: can just concatenate default_val and value
        spec.children = override_flags(arg_report(reptype,source,[default_val value skip_noreps]),reflag{:});
    end
end
