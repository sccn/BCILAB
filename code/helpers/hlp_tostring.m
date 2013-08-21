function str = hlp_tostring(v,stringcutoff,prec)
% Get an human-readable string representation of a data structure.
% String = hlp_tostring(Data,StringCutoff)
%
% The resulting string representations are usually executable, but there are corner cases (e.g.,
% certain anonymous function handles and large data sets), which are not supported. For
% general-purpose serialization, see hlp_serialize/hlp_deserialize.
%
% In:
%   Data : a data structure
%
%   StringCutoff : optional maximum string length for atomic fields (default: 0=off)
%
%   Precision : maximum significant digits (default: 15)
%
% Out:
%   String : string form of the data structure
%
% Notes:
%   hlp_tostring has builtin support for displaying expression data structures.
%
% Examples:
%   % get a string representation of a data structure
%   hlp_tostring({'test',[1 2 3], struct('field','value')})
%
% See also:
%   hlp_serialize
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15
%
%                                adapted from serialize.m
%                                (C) 2006 Joger Hansegord (jogerh@ifi.uio.no)

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

if nargin < 2
    stringcutoff = 0; end
if nargin < 3
    prec = 15; end

str = serializevalue(v);


    function val = serializevalue(v)
        % Main hub for serializing values
        if isnumeric(v) || islogical(v)
            val = serializematrix(v);
        elseif ischar(v)
            val = serializestring(v);
        elseif isa(v,'function_handle')
            val = serializefunction(v);
        elseif is_impure_expression(v)
            val = serializevalue(v.tracking.expression);
        elseif has_canonical_representation(v)
            val = serializeexpression(v);
        elseif is_dataset(v)
            val = serializedataset(v);
        elseif isstruct(v)
            val = serializestruct(v);
        elseif iscell(v)
            val = serializecell(v);
        elseif isobject(v)
            val = serializeobject(v);
        else
            try
                val = serializeobject(v);
            catch
                error('Unhandled type %s', class(v));
            end
        end
    end

    function val = serializestring(v)
        % Serialize a string
        if any(v == '''')
            val = ['''' strrep(v,'''','''''') ''''];
            try
                if ~isequal(eval(val),v)
                    val = ['char(' serializevalue(uint8(v)) ')']; end
            catch
                val = ['char(' serializevalue(uint8(v)) ')'];
            end
        else
            val = ['''' v ''''];
        end
        val = trim_value(val,'''');
    end

    function val = serializematrix(v)
        % Serialize a matrix and apply correct class and reshape if required
        if ndims(v) < 3 %#ok<ISMAT>
            if isa(v, 'double')
                if size(v,1) == 1 && length(v) > 3 && isequal(v,v(1):v(2)-v(1):v(end))
                    % special case: colon sequence
                    if v(2)-v(1) == 1
                        val = ['[' num2str(v(1)) ':' num2str(v(end)) ']'];
                    else
                        val = ['[' num2str(v(1)) ':' num2str(v(2)-v(1)) ':' num2str(v(end)) ']'];
                    end
                elseif size(v,2) == 1 && length(v) > 3 && isequal(v',v(1):v(2)-v(1):v(end))
                    % special case: colon sequence
                    if v(2)-v(1) == 1
                        val = ['[' num2str(v(1)) ':' num2str(v(end)) ']'''];
                    else
                        val = ['[' num2str(v(1)) ':' num2str(v(2)-v(1)) ':' num2str(v(end)) ']'''];
                    end
                else
                    val = mat2str(v,prec);
                end
            else
                val = mat2str(v,prec,'class');
            end
            val = trim_value(val,']');
        else
            if isa(v, 'double')
                val = mat2str(v(:),prec);
            else
                val = mat2str(v(:),prec,'class');
            end
            val = trim_value(val,']');
            val = sprintf('reshape(%s, %s)', val, mat2str(size(v)));
        end
    end

    function val = serializecell(v)
        % Serialize a cell
        if isempty(v)
            val = '{}';
            return
        end
        cellSep = ', ';
        if isvector(v) && size(v,1) > 1
            cellSep = '; ';
        end
        
        % Serialize each value in the cell array, and pad the string with a cell
        % separator.
        vstr = cellfun(@(val) [serializevalue(val) cellSep], v, 'UniformOutput', false);
        vstr{end} = vstr{end}(1:end-2);
        
        % Concatenate the elements and add a reshape if requied
        val = [ '{' vstr{:} '}'];
        if ~isvector(v)
            val = ['reshape('  val sprintf(', %s)', mat2str(size(v)))];
        end
    end

    function val = serializeexpression(v)
        % Serialize an expression
        if numel(v) > 1
            val = ['['];
            for k = 1:numel(v)
                val = [val serializevalue(v(k)), ', ']; end
            val = [val(1:end-2) ']'];
        else
            if numel(v.parts) > 0
                val = [char(v.head) '('];
                for fieldNo = 1:numel(v.parts)
                    val = [val serializevalue(v.parts{fieldNo}), ', ']; end
                val = [val(1:end-2) ')'];
            else
                val = char(v.head);
            end
        end
    end

    function val = serializedataset(v) %#ok<INUSD>
        % Serialize a data set
        val = '<EEGLAB data set>';
    end

    function val = serializestruct(v)
        % Serialize a struct by converting the field values using struct2cell
        fieldNames   = fieldnames(v);
        fieldValues  = struct2cell(v);
        if ndims(fieldValues) > 6
            error('Structures with more than six dimensions are not supported');
        end
        val = 'struct(';
        for fieldNo = 1:numel(fieldNames)
            val = [val serializevalue( fieldNames{fieldNo}) ', '];
            val = [val serializevalue( permute(fieldValues(fieldNo, :,:,:,:,:,:), [2:ndims(fieldValues) 1]) ) ];
            val = [val ', '];
        end
        if numel(fieldNames)==0
            val = [val ')'];
        else
            val = [val(1:end-2) ')'];
        end
        if ~isvector(v)
            val = sprintf('reshape(%s, %s)', val, mat2str(size(v)));
        end
    end

    function val = serializeobject(v)
        % Serialize an object by converting to struct and add a call to the copy constructor
        val = sprintf('%s(%s)', class(v), serializevalue(struct(v)));
    end


    function val = serializefunction(v)
        % Serialize a function handle
        try
            val = ['@' char(get_function_symbol(v))];
        catch
            val = char(v);
        end
    end

    function v = trim_value(v,appendix)
        if nargin < 2
            appendix = ''; end
        % Trim a serialized value to a certain length
        if stringcutoff && length(v) > stringcutoff
            v = [v(1:stringcutoff) '...' appendix]; end
    end

end

function result___ = get_function_symbol(expression___)
% internal: some function_handle expressions have a function symbol (an @name expression), and this function obtains it
% note: we are using funny names here to bypass potential name conflicts within the eval() clause further below
if ~isa(expression___,'function_handle')
    error('the expression has no associated function symbol.'); end

string___ = char(expression___);
if string___(1) == '@'
    % we are dealing with a lambda function
    if is_symbolic_lambda(expression___)
        result___ = eval(string___(27:end-21));
    else
        error('cannot derive a function symbol from a non-symbolic lambda function.');
    end
else
    % we are dealing with a regular function handle
    result___ = expression___;
end
end

function res = is_symbolic_lambda(x)
% internal: a symbolic lambda function is one which generates expressions when invoked with arguments (this is what exp_symbol generates)
res = isa(x,'function_handle') && ~isempty(regexp(char(x),'@\(varargin\)struct\(''head'',\{.*\},''parts'',\{varargin\}\)','once'));
end
