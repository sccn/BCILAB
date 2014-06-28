function spec = expand_argsubswitch(varargin)
% Internal: expand the output of an arg_subswitch(...) declaration function into an argument specifier.
% Specifier = expand_argsubswitch(ReportType, ...)
%
% The argument declaration functions used in an arg_define declaration produce a compact and
% low-overhead representation that first needs to be expanded into a full specifier struct that
% contains all properties of the argument and can be processed by other functions. This function
% performs the expansion for the arg_subswitch() arguments.
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
    spec = hlp_microcache('arg',@expand_specifier,varargin{:});
end


% expand the arg_subswitch(...) declaration line into an argument specifier (see expand_argsubtoggle for a simpler example)
function spec = expand_specifier(reptype,names,defaults,sources,help,varargin)
    % set defaults
    if nargin < 4 || isempty(sources)
        error('BCILAB:args:no_options','The Sources argument for arg_subswitch() may not be omitted.'); end 
    if iscellstr(sources(1:2:end)) && all(cellfun(@(x)iscell(x)||isa(x,'function_handle'),sources(2:2:end)))
        sources = mat2cell(sources,1,repmat(2,length(sources)/2,1)); end % we turn the NVP form of sources into a cell array of cell arrays
    if nargin < 5
        help = ''; end
    
    % extract special options unique to arg_subswitch
    [reflag,suppress,range,case_defaults] = deal({},{},cellfun(@(c)c{1},sources,'UniformOutput',false),cell(1,length(sources)));
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'reflag','suppress'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end
    
    % handle 'reflag' option (remap from {'key3',reflag3,'key1',reflag1, ...} to {reflag1,reflag2,reflag3,...})
    if ~isempty(reflag)
        keys = reflag(1:2:end);
        values = reflag(2:2:end);
        if ~isempty(fast_setdiff(keys,range))
            error(['Some of the keys in reflag are not defined in Sources: ' hlp_tostring(fast_setdiff(keys,range))]); end
        reflag = repmat({{}},1,length(range));
        for k=1:length(keys)
            reflag{strcmp(keys{k},range)} = values{k}; end
    else
        reflag = repmat({{}},1,length(range));
    end
    
    % handle the 'suppress' option (by appending to reflag)
    if ischar(suppress)
        suppress = {suppress}; end
    for n=suppress
        reflag = cellfun(@(f)[f,{n{1},{'displayable',false}}],reflag,'UniformOutput',false); end

    % reformat the Sources
    for k=length(sources):-1:1
        if ~ischar(sources{k}{1})
            error('In arg_subswitch, each selector must be a string.'); end
        % extract special arguments (Defaults and Format), if present
        fmt = [0 Inf];
        if length(sources{k}) >= 3
            case_defaults{k} = sources{k}{3}; end
        if length(sources{k}) >= 4
            fmt = sources{k}{4}; end
        % parse the k'th Source
        sources{k} = generate_source(fmt,sources{k}{2});
    end
    
    % initialize the specification struct
    spec = arg_specifier('head','arg_subswitch', 'names',names, 'help',help, 'mapper',@map_argsubswitch, varargin{:}, ...
        'range',range, 'type','char', 'shape','row', 'alternatives',cell(1,length(sources)), 'reflag',reflag, 'defaults',{defaults}, 'sources',sources);
    
    % post-process the 'mapper' (needs to take 3 arguments)
    if nargin(spec.mapper) == 1
        spec.mapper = @(x,y,z) spec.mapper(x); end

    % prepend alternative defaults (either those that are specified, or all if we're building a rich spec)
    for n=find(~cellfun('isempty',case_defaults) | strcmp(reptype,'rich'))
        spec.defaults = [{spec.range(n) case_defaults{n}} spec.defaults]; end
end


% this function maps an argument list onto a string selection key and the cell array of 
% name-value pairs / structs to assign
function [selection,args] = map_argsubswitch(args,selectors,names)
    % perform type checking
    if ~iscell(args)
        if isstruct(args) || ischar(args)
            args = {args};
        elseif isequal(args,[])
            args = {};
        else
            error(['It is not allowed to assign anything other than a cell, a struct, or a (selector) string to an arg_subswitch argument (here:' names{1} ')']); 
        end
    end
    
    % perform mapping
    if isempty(args)
        selection = selectors{1};
    elseif isfield(args{1},'arg_selection')
        selection = args{1}.arg_selection;
    elseif any(strcmp(args{1},selectors))
        [selection,args] = deal(args{1},args(2:end));
    else
        pos = find(strcmp('arg_selection',args(1:end-1)),1,'last');
        [selection,args] = deal(args{pos+1},args([1:pos-1 pos+2:end]));
    end
    
    % If this error is triggered, an value was passed for an argument which has a flexible structure (chosen out of a set of possibilities), but the possibility
    % which was chosen according to the passed value does not match any of the specified ones. For a value that is a cell array of arguments, the choice is 
    % made based on the first element in the cell. For a value that is a structure of arguments, the choice is made based on the 'arg_selection' field.
    % The error is usually resolved by reviewing the argument specification of the offending function carefully, and comparing the passed value to the Alternatives
    % declared in the arg_subswitch() clause in which the offending argument is declared.
    if isempty(selection)
        error(['The chosen selector argument (empty) does not match any of the possible options (' sprintf('%s, ',selectors{1:end-1}) selectors{end} ') in the function argument ' names{1} '.']);
    elseif ~any(strcmpi(selection,selectors))
        error(['The chosen selector argument (' selection ') does not match any of the possible options (' sprintf('%s, ',selectors{1:end-1}) selectors{end} ') in the function argument ' names{1} '.']); 
    end
end
