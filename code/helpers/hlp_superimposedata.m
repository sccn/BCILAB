function res = hlp_superimposedata(varargin)
% Merge multiple partially populated data structures into one fully populated one.
% Result = hlp_superimposedata(Data1, Data2, Data3, ...)
%
% The function is applicable when you have cell arrays or structs/struct arrays with non-overlapping
% patterns of non-empty entries, where all entries should be merged into a single data structure
% which retains their original positions. If entries exist in multiple data structures at the same
% location, entries of later items will be ignored (i.e. earlier data structures take precedence).
%
% In:
%   DataK : a data structure that should be super-imposed with the others to form a single data
%           structure
%
% Out:
%   Result : the resulting data structure
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-19

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

% first, compactify the data by removing the empty items
compact = varargin(~cellfun('isempty',varargin));
if isempty(compact)
    res = [];
else
    % start with the last data structure, then merge the remaining data structures into it (in reverse
    % order as this avoids having to grow arrays incrementally in typical cases)
    res = compact{end};
    for k=length(compact)-1:-1:1
        res = merge(res,compact{k}); end
end

% merge data structures A and B
function A = merge(A,B)

if iscell(A) && iscell(B)
    % make sure that both have the same number of dimensions
    if ndims(A) > ndims(B)
        B = grow_cell(B,size(A));
    elseif ndims(A) < ndims(B)
        A = grow_cell(A,size(B));
    end
    
    % make sure that both have the same size
    if all(size(B)==size(A))
        % we're fine
    elseif all(size(B)>=size(A))
        % A is a minor of B: grow A
        A = grow_cell(A,size(B));
    elseif all(size(A)>=size(B))
        % B is a minor of A: grow B
        B = grow_cell(B,size(A));
    else
        % A and B have mixed sizes... grow both as necessary
        M = max(size(A),size(B));
        A = grow_cell(A,M);
        B = grow_cell(B,M);
    end
    
    % find all non-empty elements in B
    idx = find(~cellfun(@(x)isequal(x,[]),B));
    
    if ~isempty(idx)
        % check if any of these is occupied in A
        clean = cellfun('isempty',A(idx));
        if ~all(clean)
            % merge all conflicting items recursively
            conflicts = idx(~clean);
            for k=conflicts(:)'
                A{k} = merge(A{k},B{k}); end
            % and transfer the rest
            if any(clean)
                A(idx(clean)) = B(idx(clean)); end
        else
            % transfer all to A
            A(idx) = B(idx);
        end
    end
elseif isstruct(A) && isstruct(B)
    % first make sure that both have the same fields
    fnA = fieldnames(A);
    fnB = fieldnames(B);
    if isequal(fnA,fnB)
        % we're fine
    elseif isequal(sort(fnA),sort(fnB))
        % order doesn't match -- impose A's order on B
        B = orderfields(B,fnA);
    elseif isempty(setdiff(fnA,fnB))
        % B has a superset of A's fields: add the remaining fields to A, and order them according to B
        remaining = setdiff(fnB,fnA);
        for fn = remaining'
            A(1).(fn{1}) = []; end
        A = orderfields(A,fnB);
    elseif isempty(setdiff(fnB,fnA))
        % A has a superset of B's fields: add the remaining fields to B, and order them according to A
        remaining = setdiff(fnA,fnB);
        for fn = remaining'
            B(1).(fn{1}) = []; end
        B = orderfields(B,fnA);
    else
        % A and B have incommensurable fields; add B's fields to A's fields, add A's fields to B's
        % and order according to A's fields
        remainingB = setdiff(fnB,fnA);
        for fn = remainingB'
            A(1).(fn{1}) = []; end
        remainingA = setdiff(fnA,fnB);
        for fn = remainingA'
            B(1).(fn{1}) = []; end
        B = orderfields(B,A);
    end
    
    % that being established, convert them to cell arrays, merge their cell arrays, and convert back to structs
    merged = merge(struct2cell(A),struct2cell(B));
    A = cell2struct(merged,fieldnames(A),1);
    
elseif isstruct(A) && ~isstruct(B)
    if ~isempty(B)
        error('One of the sub-items is a struct, and the other one is of a non-struct type.');
    else
        % we retain A
    end
elseif isstruct(B) && ~isstruct(A)
    if ~isempty(A)
        error('One of the sub-items is a struct, and the other one is of a non-struct type.');
    else
        % we retain B
        A = B;
    end
elseif iscell(A) && ~iscell(B)
    if ~isempty(B)
        error('One of the sub-items is a cell array, and the other one is of a non-cell type.');
    else
        % we retain A
    end
elseif iscell(B) && ~iscell(A)
    if ~isempty(A)
        error('One of the sub-items is a cell array, and the other one is of a non-cell type.');
    else
        % we retain B
        A = B;
    end
elseif isempty(A) && ~isempty(B)
    % we retain B
    A = B;
elseif isempty(B) && ~isempty(A)
    % we retain A
elseif ~isequal_weak(A,B)
    % we retain A and warn about dropping B
    warn_once('Two non-empty (and non-identical) sub-elements occupied the same index; one was dropped. This warning will only be displayed once.');
end

% grow a cell array to accomodate a particular index
% (assuming that this index is not contained in the cell array yet)
function C = grow_cell(C,idx)
tmp = sprintf('%i,',idx);
eval(['C{' tmp(1:end-1) '} = [];']);
