function [closest,rest] = utl_collection_closest(varargin)
% Find the best-matching data set(s) in a collection for a given identifier record
%
% In:
%   Collection : data set collection
%
%   Identifier : an identifier record, containing meta-data; either a struct with fields or a cell
%                array with name-value pairs
%
%   ScopeOrdering : Odering of matching-relevant properties from largest difference to smallest
%                   (default: {'group','subject','day','montage','session','recording','block'}).
%
%   TreatAsFlat : Cell array of property names to treat as flat. Normally, the closest session to 
%                 session 1 would be session 2; however, if other sessions shall have the same
%                 distance, then ''session'' should be added to this list. This will allow to sort,
%                 for instance, the data from all other sessions into the "best-matching" set, since
%                 they will all have the same minimum distance. (default: {'group','subject',
%                 'montage','block'})
%
% Out:
%   Result : the set in the collection whose meta-data best matches that of the identifier,
%            or multiple sets if their distance to the identifier is identical
%
% See also:
%   bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-29
dp;

args = arg_define(0:2, varargin, ...
    arg({'collection','Collection'},mandatory,[],'A dataset collection. Cell array of structs with meta-data.','type','expression'), ...
    arg({'identifier','Identifier'},mandatory,[],'An identifier record, containing meta-data; either a struct with fields or a cell array with name-value pairs.','type','expression'), ...
    arg({'scope_order','ScopeOrdering'},{'group','subject','day','montage','session','recording','block'},[],'Dataset identifiers ordered from coarsest to finest. This both defines what identifiers are known, and their hierarchical ordering.'), ...    
    arg({'treat_as_flat','TreatAsFlat'},{'group','subject','montage','block'},[],'Cell array of property names to treat as flat. Normally, the closest session to session 1 would be session 2; however, if other sessions shall have the same distance, then ''session'' should be added to this. This will allow to sort, for instance, the data from all other sessions into the "best-matching" set.','type','expression'));
    
[collection,identifier] = deal(args.collection, args.identifier);


if ~iscell(collection) || ~all(cellfun('isclass',collection,'struct'))
    error('The given Collection argument must be a cell array of structs.'); end

if ~isempty(identifier)
    if iscell(identifier) && iscellstr(identifier(1:2:end))
        identifier = hlp_varargin2struct(identifier); end
    if ~isscalar(identifier) || ~isstruct(identifier)
        error('The identifier should be a 1x1 struct'); end
    
    % our default ordering hierarchy of dis-similarity (first is largest granularity, last is finest)
    args.scope_order = {'subject','day','montage','session','recording','block'};
    
    
    % first tag all collection items with a tracking index
    for k=1:length(collection)
        collection{k}.tracking_index = k; end
    original_collection = collection;
    
    
    % first go through each scale order and restrict the recordings appropriately
    for s = 1:length(args.scope_order)
        id = args.scope_order{s};
        % is the order present in the identifier?
        if isfield(identifier,id)
            present = find(cellfun(@(x)isfield(x,id),collection));
            % is it present in the collection?
            if ~isempty(present)
                values = cellfun(@(x)x.(id),collection(present),'UniformOutput',false);
                if iscellstr(values)
                    % if it's a string value, is there a perfect match?
                    perfect_matches = strcmp(values,identifier.(id));
                    if any(perfect_matches)
                        % then restrict the collection to that
                        collection = collection(present(perfect_matches));
                    else
                        fprintf('Note: no element in the collection matches the value "%s" of the identifier''s "%s" field.\n',identifier.(id),id);
                    end
                elseif all(cellfun(@isnumeric,values))
                    % if it's a numeric value, retain those that are closest to the identifier
                    if any(strcmp(id, args.treat_as_flat))
                        distance = double(cellfun(@(x)x,values) ~= identifier.(id));
                    else
                        distance = abs(cellfun(@(x)x,values) - identifier.(id));
                    end
                    retain = distance == min(distance);
                    collection = collection(present(retain));
                else
                    % otherwise use isequal_weak for perfect matches
                    perfect_matches = cellfun(@(x)isequal_weak(x,identifier.(id)),values);
                    if any(perfect_matches)
                        % then restrict the collection to that
                        collection = collection(present(perfect_matches));
                    else
                        fprintf('Note: no element in the collection matches the value "%s" of the identifier''s "%s" field.\n',hlp_tostring(identifier.(id)),id);
                    end
                end
            end
        end
    end
    
    % find remaining fields and values of the identifier
    rest_fields = setdiff(fieldnames(identifier),args.scope_order);
    
    if ~isempty(rest_fields)
        rest_values = cellfun(@(id)identifier.(id),rest_fields,'UniformOutput',false);
        
        % for each remaining non-numeric field in the identifier...
        numeric = cellfun(@isnumeric,rest_values);
        for f=find(~numeric)
            id = rest_fields{f};
            present = find(cellfun(@(x)isfield(x,id),collection));
            % is it present in the collection?
            if ~isempty(present)
                values = cellfun(@(x)x.(id),collection(present),'UniformOutput',false);
                % find perfect matches
                if iscellstr(values)
                    perfect_matches = strcmp(values,identifier.(id));
                else
                    perfect_matches = cellfun(@(x)isequal_weak(x,identifier.(id)),values);
                end
                if any(perfect_matches)
                    % then restrict the collection to that
                    collection = collection(present(perfect_matches));
                else
                    fprintf('Note: no element in the collection matches the value "%s" of the identifier''s "%s" field.\n',hlp_tostring(identifier.(id)),id);
                end
            end
        end
        
        % get rid of those numeric fields that are not present in the collection
        % and get rid of collection fields that are lacking a numeric field that other items have
        retain_numeric = [];
        for f=find(numeric)
            id = rest_fields{f};
            present = find(cellfun(@(x)isfield(x,id),collection));
            if ~isempty(present)
                % retain the numeric field
                retain_numeric(end+1) = f; %#ok<AGROW>
                % remove collection items that are lacking this field
                collection = collection(present);
            end
        end
        
        % concatenate the remaining numeric properties of the collection items
        if ~isempty(retain_numeric)
            collection_data = [];
            identifier_data = [];
            for f=retain_numeric
                id = rest_fields{f};
                values = cellfun(@(x)x.(id)(:),collection);
                collection_data = [collection_data; values]; %#ok<AGROW>
                identifier_data = [identifier_data; identifier.(id)(:)]; %#ok<AGROW>
            end
            
            % compute euclidean distances
            distances = sum(sqr(bsxfun(@minus,collection_data,identifier_data)));
            retain = distances == min(distances);
            collection = collection(retain);
        end
    end
    
    % finish the set of closest matches
    closest = collection;
    closest_indices = cellfun(@(x)x.tracking_index,closest);
    for k=1:length(closest)
        closest{k} = rmfield(closest{k},'tracking_index'); end
    
    % finish the set of remaining sets
    not_closest = cellfun(@(x)~ismember(x.tracking_index,closest_indices),original_collection);
    rest = original_collection(not_closest);
    for k=1:length(rest)
        rest{k} = rmfield(rest{k},'tracking_index'); end
    
else
    closest = collection;
    rest = {};
end