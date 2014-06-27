function spec = expand_argsub(unused,varargin) %#ok<INUSL>
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
    spec = hlp_microcache('arg',@expand_specifier,varargin{:});
end


% expand the arg_sub(...) declaration line into an argument specifier (see expand_arg for a simpler example)
function spec = expand_specifier(names,defaults,source,help,varargin)
    % set defaults
    if nargin < 1
        names = {}; end
    if nargin < 2
        defaults = {}; end
    if nargin < 3
        source = {}; end
    if nargin < 4
        help = ''; end    
    
    % extract special options unique to arg_sub
    [fmt,suppress] = deal([0 Inf],{});
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'fmt','suppress'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end    
    
    % initialize the specification struct
    spec = arg_specifier('head','arg_sub', 'names',names, 'help',help, 'mapper',@map_argsub, varargin{:}, ...
        'value',[], 'type','char', 'shape','row', 'sources',{generate_source(fmt,source)}, 'defaults',{defaults});

    % post-process the 'mapper' (needs to take 3 arguments)
    if ~isempty(spec.mapper) && nargin(spec.mapper) == 1
        spec.mapper = @(x,y,z) spec.mapper(x); end
    
    % handle the 'suppress' option (by appending to reflag)
    if ischar(suppress)
        suppress = {suppress}; end
    for n=suppress
        spec.reflag = [spec.reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % post-process the 'reflag' option (wrapped in a cell for consistency with arg_subswitch)
    spec.reflag = {spec.reflag};    
end

% this function maps an argument list onto a selection key and the cell array of 
% name-value pairs / structs to assign; for arg_sub the key is empty
function [selection,args] = map_argsub(args,varargin)
    if ~iscell(args)
        args = {args}; end
    selection = [];
end
