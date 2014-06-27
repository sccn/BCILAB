function v = hlp_deserialize(m)
% Convert a serialized byte vector back into the corresponding MATLAB data structure.
% Data = hlp_deserialize(Bytes)
%
% In:
%   Bytes : a representation of the original data as a byte stream
%
% Out:
%   Data : some MATLAB data structure
%
%
% See also:
%   hlp_serialize
%
% Examples:
%   bytes = hlp_serialize(mydata);
%   ... e.g. transfer the 'bytes' array over the network ...
%   mydata = hlp_deserialize(bytes);
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-02
%
%                                adapted from deserialize.m
%                                (C) 2010 Tim Hutt

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

% wrap dispatcher
v = deserialize_value(uint8(m(:)),1);

end

% dispatch
function [v,pos] = deserialize_value(m,pos)
switch m(pos)
    case {0,200}
        [v,pos] = deserialize_string(m,pos);
    case 128
        [v,pos] = deserialize_struct(m,pos);
    case {33,34,35,36,37,38,39}
        [v,pos] = deserialize_cell(m,pos);
    case {1,2,3,4,5,6,7,8,9,10}
        [v,pos] = deserialize_scalar(m,pos);
    case 133
        [v,pos] = deserialize_logical(m,pos);
    case {151,152,153}
        [v,pos] = deserialize_handle(m,pos);
    case {17,18,19,20,21,22,23,24,25,26}
        [v,pos] = deserialize_numeric_simple(m,pos);
    case 130
        [v,pos] = deserialize_sparse(m,pos);
    case 131
        [v,pos] = deserialize_complex(m,pos);
    case 132
        [v,pos] = deserialize_char(m,pos);
    case 134
        [v,pos] = deserialize_object(m,pos);
    otherwise
        error('Unknown class');
end
end

% individual scalar
function [v,pos] = deserialize_scalar(m,pos)
classes = {'double','single','int8','uint8','int16','uint16','int32','uint32','int64','uint64'};
sizes = [8,4,1,1,2,2,4,4,8,8];
sz = sizes(m(pos));
% Data.
v = typecast(m(pos+1:pos+sz),classes{m(pos)});
pos = pos + 1 + sz;
end

% standard string
function [v,pos] = deserialize_string(m,pos)
if m(pos) == 0
    % horizontal string: tag
    pos = pos + 1;
    % length (uint32)
    nbytes = double(typecast(m(pos:pos+3),'uint32'));
    pos = pos + 4;
    % data (chars)
    v = char(m(pos:pos+nbytes-1))';
    pos = pos + nbytes;
else
    % proper empty string: tag
    [v,pos] = deal('',pos+1);
end
end

% general char array
function [v,pos] = deserialize_char(m,pos)
pos = pos + 1;
% Number of dims
ndms = double(m(pos));
pos = pos + 1;
% Dimensions
dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
pos = pos + ndms*4;
nbytes = prod(dms);
% Data.
v = char(m(pos:pos+nbytes-1));
pos = pos + nbytes;
v = reshape(v,[dms 1 1]);
end

% general logical array
function [v,pos] = deserialize_logical(m,pos)
pos = pos + 1;
% Number of dims
ndms = double(m(pos));
pos = pos + 1;
% Dimensions
dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
pos = pos + ndms*4;
nbytes = prod(dms);
% Data.
v = logical(m(pos:pos+nbytes-1));
pos = pos + nbytes;
v = reshape(v,[dms 1 1]);
end

% simple numerical matrix
function [v,pos] = deserialize_numeric_simple(m,pos)
classes = {'double','single','int8','uint8','int16','uint16','int32','uint32','int64','uint64'};
sizes = [8,4,1,1,2,2,4,4,8,8];
cls = classes{m(pos)-16};
sz = sizes(m(pos)-16);
pos = pos + 1;
% Number of dims
ndms = double(m(pos));
pos = pos + 1;
% Dimensions
dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
pos = pos + ndms*4;
nbytes = prod(dms) * sz;
% Data.
v = typecast(m(pos:pos+nbytes-1),cls);
pos = pos + nbytes;
v = reshape(v,[dms 1 1]);
end

% complex matrix
function [v,pos] = deserialize_complex(m,pos)
pos = pos + 1;
[re,pos] = deserialize_numeric_simple(m,pos);
[im,pos] = deserialize_numeric_simple(m,pos);
v = complex(re,im);
end

% sparse matrix
function [v,pos] = deserialize_sparse(m,pos)
pos = pos + 1;
% matrix dims
u = double(typecast(m(pos:pos+7),'uint64'));
pos = pos + 8;
v = double(typecast(m(pos:pos+7),'uint64'));
pos = pos + 8;
% index vectors
[i,pos] = deserialize_numeric_simple(m,pos);
[j,pos] = deserialize_numeric_simple(m,pos);
if m(pos)
    % real
    pos = pos+1;
    [s,pos] = deserialize_numeric_simple(m,pos);
else
    % complex
    pos = pos+1;
    [re,pos] = deserialize_numeric_simple(m,pos);
    [im,pos] = deserialize_numeric_simple(m,pos);
    s = complex(re,im);
end
v = sparse(i,j,s,u,v);
end

% struct array
function [v,pos] = deserialize_struct(m,pos)
pos = pos + 1;
% Number of field names.
nfields = double(typecast(m(pos:pos+3),'uint32'));
pos = pos + 4;
% Field name lengths
fnLengths = double(typecast(m(pos:pos+nfields*4-1),'uint32'));
pos = pos + nfields*4;
% Field name char data
fnChars = char(m(pos:pos+sum(fnLengths)-1)).';
pos = pos + length(fnChars);
% Number of dims
ndms = double(typecast(m(pos:pos+3),'uint32'));
pos = pos + 4;
% Dimensions
dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
pos = pos + ndms*4;
% Field names.
fieldNames = cell(length(fnLengths),1);
splits = [0; cumsum(double(fnLengths))];
for k=1:length(splits)-1
    fieldNames{k} = fnChars(splits(k)+1:splits(k+1)); end
% Content.
v = reshape(struct(),[dms 1 1]);
if m(pos)
    % using struct2cell
    pos = pos + 1;
    [contents,pos] = deserialize_cell(m,pos);
    v = cell2struct(contents,fieldNames,1);
else
    % using per-field cell arrays
    pos = pos + 1;
    for ff = 1:nfields
        [contents,pos] = deserialize_cell(m,pos);
        [v.(fieldNames{ff})] = deal(contents{:});
    end
end
end

% cell array
function [v,pos] = deserialize_cell(m,pos)
kind = m(pos);
pos = pos + 1;
switch kind
    case 33 % arbitrary/heterogenous cell array
        % Number of dims
        ndms = double(m(pos));
        pos = pos + 1;
        % Dimensions
        dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
        pos = pos + ndms*4;
        % Contents
        v = cell([dms,1,1]);
        for ii = 1:numel(v)
            [v{ii},pos] = deserialize_value(m,pos); end
    case 34 % cell scalars
        [content,pos] = deserialize_value(m,pos);
        v = cell(size(content));
        for k=1:numel(v)
            v{k} = content(k); end
    case 35 % mixed-real cell scalars
        [content,pos] = deserialize_value(m,pos);
        v = cell(size(content));
        for k=1:numel(v)
            v{k} = content(k); end
        [reality,pos] = deserialize_value(m,pos);
        v(reality) = real(v(reality));
    case 36 % cell array with horizontal or empty strings
        [chars,pos] = deserialize_string(m,pos);
        [lengths,pos] = deserialize_numeric_simple(m,pos);
        [empty,pos] = deserialize_logical(m,pos);
        v = cell(size(lengths));
        splits = [0 cumsum(double(lengths(:)))'];
        for k=1:length(lengths)
            v{k} = chars(splits(k)+1:splits(k+1)); end
        [v{empty}] = deal('');
    case 37 % empty,known type
        tag = m(pos);
        pos = pos + 1;
        switch tag
            case 1   % double - []
                prot = [];
            case 33  % cell - {}
                prot = {};
            case 128 % struct - struct()
                prot = struct([]);
            otherwise
                error('Unsupported type tag.');
        end
        % Number of dims
        ndms = double(m(pos));
        pos = pos + 1;
        % Dimensions
        dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
        pos = pos + ndms*4;
        % Create content
        v = repmat({prot},dms);
    case 38 % empty, prototype available
        % Prototype.
        [prot,pos] = deserialize_value(m,pos);
        % Number of dims
        ndms = double(m(pos));
        pos = pos + 1;
        % Dimensions
        dms = double(typecast(m(pos:pos+ndms*4-1),'uint32')');
        pos = pos + ndms*4;
        % Create content
        v = repmat({prot},dms);
    case 39 % boolean flags
        [content,pos] = deserialize_logical(m,pos);
        v = cell(size(content));
        for k=1:numel(v)
            v{k} = content(k); end
    otherwise
        error('Unsupported cell array type.');
end
end

% object
function [v,pos] = deserialize_object(m,pos)
pos = pos + 1;
% Get class name.
[cls,pos] = deserialize_string(m,pos);
% Get contents
[conts,pos] = deserialize_value(m,pos); 
% construct object
try
    % try to use the loadobj function
    v = eval([cls '.loadobj(conts)']);
catch
    try
        % pass the struct directly to the constructor
        v = eval([cls '(conts)']);
    catch
        try
            % try to set the fields manually
            v = feval(cls);
            for fn=fieldnames(conts)'
                try
                    set(v,fn{1},conts.(fn{1})); 
                catch
                    % Note: if this happens, your deserialized object might not be fully identical
                    % to the original (if you are lucky, it didn't matter, through). Consider 
                    % relaxing the access rights to this property or add support for loadobj from
                    % a struct.
                    warn_once('hlp_deserialize:restricted_access','No permission to set property %s in object of type %s.',fn{1},cls);
                end
            end
        catch
            v = conts;
            v.hlp_deserialize_failed = ['could not construct class: ' cls];
        end
    end
end
end

% function handle
function [v,pos] = deserialize_handle(m,pos)
% Tag
kind = m(pos);
pos = pos + 1;
switch kind
    case 151 % simple function
        persistent db_simple; %#ok<TLEV> % database of simple functions (indexed by name)
        % Name
        [name,pos] = deserialize_string(m,pos);
        try
            % look up from table
            v = db_simple.(name);
        catch
            % otherwise generate & fill table
            v = str2func(name);
            db_simple.(name) = v;
        end
    case 152 % anonymous function
        % Function code
        [code,pos] = deserialize_string(m,pos);
        % Workspace
        [wspace,pos] = deserialize_struct(m,pos);
        % Construct
        v = restore_function(code,wspace);
    case 153 % scoped or nested function
        persistent db_nested; %#ok<TLEV> % database of nested functions (indexed by name)
        % Parents
        [parentage,pos] = deserialize_cell(m,pos);
        try
            key = sprintf('%s_',parentage{:});
            % look up from table
            v = db_nested.(key);
        catch
            % recursively look up from parents, assuming that these support the arg system
            v = parentage{end};
            for k=length(parentage)-1:-1:1
                % Note: if you get an error here, you are trying to deserialize a function handle
                % to a nested function. This is not natively supported by MATLAB and can only be made
                % to work if your function's parent implements some mechanism to return such a handle.
                % The below call assumes that your function uses the BCILAB arg system to do this.
                try
                    next_v = arg_report('handle',v,parentage(k));
                catch
                    warn_once('hlp_deserialize:lookup_failed',['Could not look report properties of scoped/nested function handle "' parentage{k} '" from enclosing function "' char(v) '".']);
                    v = @error_deserializing_function;
                    return
                end
                if isempty(next_v{1})
                    warn_once('hlp_deserialize:lookup_failed',['Could not look up scoped/nested function handle "' parentage{k} '" from enclosing function "' char(v) '".']);
                end
                v = next_v{1};
            end
            if ~isempty(v)
                db_nested.(key) = v; 
            end
        end
end
end

% helper for deserialize_handle
function f = restore_function(decl__,workspace__)
% create workspace
for fn__=fieldnames(workspace__)'
    % we use underscore names here to not run into conflicts with names defined in the workspace
    eval([fn__{1} ' = workspace__.(fn__{1}) ;']); 
end
clear workspace__ fn__;
% evaluate declaration
f = eval(decl__);
end

% emit a specific warning only once (per MATLAB session)
function warn_once(varargin)
persistent displayed_warnings;
% determine the message content
if length(varargin) > 1 && any(varargin{1}==':') && ~any(varargin{1}==' ') && ischar(varargin{2})
    message_content = [varargin{1} sprintf(varargin{2:end})];
else
    message_content = sprintf(varargin{1:end});
end
% generate a hash of of the message content
str = java.lang.String(message_content);
message_id = sprintf('x%.0f',str.hashCode()+2^31);
% and check if it had been displayed before
if ~isfield(displayed_warnings,message_id)
    % emit the warning
    warning(varargin{:});
    % remember to not display the warning again
    displayed_warnings.(message_id) = true;
end
end
