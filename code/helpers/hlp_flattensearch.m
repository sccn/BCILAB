function x = hlp_flattensearch(x,form)
% Flatten search() clauses in a nested data structure into a flat search() clause. 
% Result = hlp_flattensearch(Expression, Output-Form)
%
% Internal tool used by utl_gridsearch to enable the specification of search parameters using
% search() clauses.
%
% In:
%   Expression  : some data structure, usually an argument to utl_gridsearch, may or may not contain
%                 nested search clauses.
%
%   Output-Form : form of the output (default: 'search')
%                  * 'search': the output shall be a flattened search clause (or a plain value if no 
%                              search)
%                  * 'cell': the output shall be a cell array of elements to search over
%
% Out:
%   Result      : a flattened search clause (or plain value), or a cell array of search possibilities.
% 
% See also:
%   search, utl_gridsearch
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-29

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


x = flatten(x);

if ~exist('form','var') || isempty(form) || strcmp(form,'search')
    % turn from cell format into search format
    if isscalar(x)
        % just one option: return the plain value
        x = x{1};
    else
        % multiple options: return a search expression
        x = struct('head',{@search},'parts',{x});
    end
elseif ~strcmp(form,'cell')
    error(['Unsupported output form: ' form]);
end


% recursively factor search expressions out of a data structure, to give an overall search
% expression (output format is cell array of search options)
function parts = flatten(x)
if isstruct(x)
    if isfield(x,{'data','srate','chanlocs','event','epoch'})
        % data set? do not descend further
        parts = {x};
    elseif all(isfield(x,{'head','parts'})) && numel(x)==1 && strcmp(char(x.head),'search')
        % search expression: flatten any nested searches...
        parts = cellfun(@flatten,x.parts,'UniformOutput',false);
        % ... and splice their parts in
        parts = [parts{:}];
    else
        % generic structure: create a cartesian product over field-wise searches
        if isscalar(x)
            parts = {x};
            % flatten per-field contents
            fields = cellfun(@flatten,struct2cell(x),'UniformOutput',false);
            lengths = cellfun('length',fields);
            % was any one a search?
            if any(lengths>1)
                fnames = fieldnames(x);
                % for each field that is a search...
                for k=find(lengths>1)'
                    % replicate all parts once for each search item in the current field
                    partnum = length(parts);
                    parts = repmat(parts,1,lengths(k));
                    % and fill each item into the appropriate place
                    for j=1:length(parts)
                        parts{j}.(fnames{k}) = fields{k}{ceil(j/partnum)}; end
                end
            end
        elseif ~isempty(x)
            % struct array (either with nested searches or a concatenation of search() expressions):
            % handle as a cell array of structs
            parts = flatten(arrayfun(@(s){s},x));
            % got a search?
            if ~isscalar(parts)
                % re-concatenate the cell contents of each part of the search expression into 
                % struct arrays
                for i=1:length(parts)
                    parts{i} = reshape([parts{i}{:}],size(parts{i})); end
            else
                parts = {x};
            end
        else
            parts = {x};
        end
    end
elseif iscell(x)
    % cell array: create a cartesian product over cell-wise searches
    parts = {x};
    x = cellfun(@flatten,x,'UniformOutput',false);
    % for each cell that is a search...
    for c=find(cellfun('length',x(:)')>1)
        % replicate all parts once for each search item in the current cell
        partnum = length(parts);
        parts = repmat(parts,1,length(x{c}));
        % and fill in the new item in the appropriate place
        for j=1:length(parts)
            parts{j}{c} = x{c}{ceil(j/partnum)}; end
    end
else
    % anything else: wrap
    parts = {x};
end
