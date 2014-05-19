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
    spec = hlp_microcache('arg',@expand_specifier,varargin{:});
end


% expand the arg_subtoggle(...) declaration line into an argument specifier (see expand_argsub for a simpler example)
function spec = expand_specifier(reptype,names,defaults,source,help,varargin)
    % set defaults
    if nargin < 4 || isempty(source)
        error('BCILAB:args:no_options','The Source argument for arg_subtoggle() may not be omitted.'); end 
    if nargin < 5
        help = ''; end
    
    % extract special options unique to arg_subtoggle
    [fmt,suppress,alternative_defaults] = deal([0 Inf],{},{});
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'fmt','suppress','alternative_defaults'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end
    
    % initialize the specification struct
    spec = arg_specifier('head','arg_subtoggle', 'names',names, 'help',help, 'mapper',@map_argsubtoggle, varargin{:}, ...
        'type','logical', 'shape','scalar', 'alternatives',cell(1,2), 'sources',{[],generate_source(fmt,source)});

    % handle the 'suppress' option (suppress appends to reflag)
    if ischar(suppress)
        suppress = {suppress}; end
    for n=suppress
        spec.reflag = [spec.reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % post-process the 'reflag' option (wrapped in a cell for consistency with arg_subswitch)
    spec.reflag = {[],spec.reflag};
    
    % post-process the 'mapper' option (needs to take 3 arguments)
    if nargin(spec.mapper) == 1
        spec.mapper = @(x,y,z) spec.mapper(x); end
    
    % set up the sequence of defaults to apply
    if strcmp(reptype,'rich') || ~isempty(alternative_defaults)
        % in rich mode, we need to also assign the 'on' state ({})
        % if alternative_defaults is non-empty, we need to assign that, too
        % both are applied before the actual defaults are to be applied
        spec.defaults = {alternative_defaults,defaults};
    else
        spec.defaults = {defaults};
    end
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
    elseif isequal(args,false)
        selected = false; % FIXME: this is deprecated
        args = {};
    else
        selected = true;
        args = {args};
    end
end
