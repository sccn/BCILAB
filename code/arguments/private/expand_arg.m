function spec = expand_arg(unused,varargin) %#ok<INUSL>
% Internal: expand the output of an arg(...) declaration function into an argument specifier.
% Specifier = expand_arg(ReportType, ...)
%
% The argument declaration functions used in an arg_define declaration produce a compact and
% low-overhead representation that first needs to be expanded into a full specifier struct that
% contains all properties of the argument and can be processed by other functions. This function
% performs the expansion for the atomic arguments, including arg(), arg_nogui(), and arg_norep().
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


% expand the arg(...) declaration line into an argument specifier
function spec = expand_specifier(names,default,range,help,varargin)
    % set defaults
    if nargin < 1
        names = {}; end
    if nargin < 2
        default = []; end
    if nargin < 3
        range = []; end
    if nargin < 4
        help = ''; end
    
    % initialize specification struct
    spec = arg_specifier('head','arg', 'names',names, 'range',range, 'help',help, 'to_double',[], varargin{:}, 'mapper',[]);
    
    % parse the type
    if isempty(spec.type)
        % if the type is unspecified it is deduced 
        if iscellstr(spec.range) && ((iscellstr(default) && isempty(fast_setdiff(default,spec.range))) || (isa(default,'logical') && isscalar(default)))
            % if both the value and the range are cell-string arrays and the value is a subset of
            % the range, the type is 'logical'; alternatively the default may be a logical scalar,
            % which stands for all-true/all-false
            spec.type = 'logical';
        elseif ~isempty(default) && ~(isequal(default,'__arg_mandatory__') || isequal(default,'__arg_unassigned__'))
            % if the value is non-empty, the type is deduced from the value
            spec.type = deduce_type(default);
        elseif isnumeric(spec.range) && ~isempty(spec.range)
            % if the range is numeric and non-empty, the type is deduced from the range
            spec.type = deduce_type(spec.range);
        elseif iscell(spec.range) && isscalar(unique(cellfun(@deduce_type,spec.range,'UniformOutput',false)))
            % if the range is a cell array of uniformly-typed values, the type is the unique
            % type of the cell entries
            spec.type = deduce_type(spec.range{1});
        elseif isempty(default) && isempty(spec.range) && ~isa(default,'double')
            % if both default and range are empty, and the default is not [], we deduce the type from the default
            spec.type = deduce_type(default);
        else
            % if no other rule applies, the type is expression
            spec.type = 'expression';
        end
    end    
    
    % parse the shape    
    if isempty(spec.shape)
        % if the shape is unspecified, it is deduced
        if ~isempty(default) && ~(isequal(default,'__arg_mandatory__') || isequal(default,'__arg_unassigned__'))
            % if the value is not empty, the shape is deduced from the value
            if (isequal(default,true) || isequal(default,false) || iscellstr(default)) && iscellstr(spec.range)
                % except in the special case where the default value is true, false or a cell array of strings, 
                % while the range is a cell array of strings, the shape is set to row
                spec.shape = 'row';
            else
                spec.shape = deduce_shape(default);
            end
        elseif ischar(default) && isempty(default)
            % if the value is an empty char array, the shape is by default 'row'
            spec.shape = 'row';
        elseif ~isempty(spec.range) && isnumeric(spec.range)
            % if a non-empty numeric range is given, the shape is assumed to be a scalar
            spec.shape = 'scalar';
        elseif iscell(spec.range) && isscalar(unique(cellfun(@deduce_shape,spec.range,'UniformOutput',false)))
            % if a cell-array range with uniform shapes is given, the shape is by default the
            % unique shape of the cell entries
            spec.shape = deduce_shape(spec.range{1});
        else
            % if no other rule applies, the shape is matrix
            spec.shape = 'matrix';
        end
    end
    
    % parse the default value
    if isequal(default,[]) && iscellstr(spec.range)
        % if the value is [] but the range is a cell-string array (as in a multi-option argument),
        % the value is set to the first option
        default = spec.range{1};
    elseif islogical(default) && isscalar(default) && iscellstr(spec.range)    
        % if the value is true/false, and the range is a cell-string array (as in a set of
        % options), false maps to the empty set and true maps to the full set
        default = quickif(default,spec.range,{});
    elseif isequal(default,[])    
        % if the value is [], it is converted to the declared type
        switch spec.type
            case {'cellstr','cell'}
                default = {};
            case {'char'}
                default = '';
            case {'denserealdouble','densecomplexdouble'}
                default = full(double(default));
            case {'sparserealdouble','sparsecomplexdouble'}
                default = sparse(double(default));
            case {'denserealsingle','densecomplexsingle'}
                default = full(single(default));
            case {'sparserealsingle','sparsecomplexsingle'}
                default = sparse(single(default));
            case {'int8','uint8','int16','uint16','int32','uint32','int64','logical'}
                default = cast(default, spec.type);
        end
    end
    
    % if to_double is still undecided, enable it only if the input type is an integer
    if isempty(spec.to_double)
        spec.to_double = any(strcmp(spec.type,{'int8','uint8','int16','uint16','int32','uint32','int64'})); end

    % set up the sequence of defaults
    spec.defaults = {default};
end
