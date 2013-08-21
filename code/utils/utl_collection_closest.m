function [closest,rest] = utl_collection_closest(collection,identifier,varargin)
% find the best-matching data set(s) in a collection for a given identifier record
%
% In:
%   Collection : data set collection
%
%   Identifier : an identifier record, containing meta-data
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


if ~isempty(identifier)
    if ~isscalar(identifier) || ~isstruct(identifier)
        error('The identifier should be a 1x1 struct'); end
    
    % our default ordering hierarchy of dis-similarity (first is largest granularity, last is finest)
    scale_orders = {'subject','day','montage','session','recording'};
    
    
    % first tag all collection items with a tracking index
    for k=1:length(collection)
        collection{k}.tracking_index = k; end
    original_collection = collection;
    
    
    % first go through each scale order and restrict the recordings appropriately
    for s = 1:length(scale_orders)
        id = scale_orders{s};
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
                    distance = abs(cellfun(@(x)x,values) - identifier.(id));
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
    rest_fields = setdiff(fieldnames(identifier),scale_orders);
    
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
                retain_numeric(end+1) = f;
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
                collection_data = [collection_data; values];
                identifier_data = [identifier_data; identifier.(id)(:)];
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