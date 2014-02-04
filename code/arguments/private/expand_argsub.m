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
    [fmt,reflag,suppress] = deal([0 Inf],{},{});
    
    % extract special options
    for k=length(varargin)-1:-2:1
        if any(strcmp(varargin{k},{'fmt','reflag','suppress'}))
            eval([varargin{k} ' = varargin{k+1}; varargin([k k+1]) = [];']); end
    end    
    
    % handle the 'suppress' option (by appending to reflag)
    if ischar(suppress)
        suppress = {suppress}; end
    for n=suppress
        reflag = [reflag {n{1},{'displayable',false}}]; end %#ok<AGROW>
    
    % initialize the specification struct
    spec = arg_specifier('head',@arg_sub, 'names',names, 'help',help,varargin{:}, ...
        'value',[], 'type','char', 'shape','row');
    
    % parse the Source
    source = generate_source(fmt,source);
    
    % set up the assigner
    spec.assigner = @(spec,value) assign_argsub(spec,value,reptype,source,reflag);
    
    % assign the defaults
    spec = spec.assigner(spec,defaults);
end


% function to perform the value assignment
function spec = assign_argsub(spec,invalue,reptype,source,reflag)
    % make sure that invalue is a cell array
    if ~iscell(invalue)
        invalue = {invalue}; end
    % parse the value into a spec struct array
    value = arg_report('parse',source,[invalue {'__arg_skip__',true}]);
    % selectively override children with the value
    spec.children = override_fields(spec.children,value);
    % override flags
    spec.children = override_flags(spec.children,reflag{:});
end
