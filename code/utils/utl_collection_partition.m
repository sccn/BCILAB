function res = utl_collection_partition(varargin)
% Partition a dataset collection according to some settings
% Result = utl_collection_partition(Collection,IndexSet,Settings)
%
% This function is a valid partitioner for utl_crossval, utl_nested_crossval, and utl_searchmodel,
% which accepts a collection (cell array) of datasets. The datasets may have additional meta-data
% that is used to determine the desired partitioning (see Settings).
%
% In:
%   Collection : a dataset collection, i.e. cell array of structs with meta-data
%
%   IndexSet : the index set to use for partitioning; this partitioner is somewhat special;
%              * if this is [], the function generates the partitioning, i.e. returns
%                a cell array of partitions, each of which is a cell array of:
%                {training collection, test collection, {'goal_identifier',identifier for test collection}}
%
%              * otherwise, we assume that this is one of the partitioned collections that we 
%                generated, and just return it unmodified
%
%   Settings : the settings for the partitioning; may have the following fields (or names, if given 
%              as a cell array of NvPs):
%
%               'restrict_to' : struct with properties to which the test sets under consideration are 
%                               restricted (default: blank / no restriction)
%               'exclude'       struct with a property signature to exclude from the test sets
%                               or a cell array of multiple signatures to exclude
%               'test_on'     : identifier that determines the items to use as test sets (default:
%                               coarsest granularity out of the known identifiers or element-by-element, 
%                               if none known)
%               'consider'    : modifier that determines how to pack and optionally prune the test 
%                               items ('each', 'all', 'last', or 'allexceptfirst')
%                               (default: 'each')
%               'per'         : if consider is not each, determines the context to which the 
%                               packing and pruning operation of 'consider' refers 
%                               (default: next coarser granularity than test_on)
%
%               'scope_order' : odering of partitioning-relevant properties from coarsest granularity 
%                               to finest (default: {'group','subject','day','montage','session','recording','block'} 
%
%                               This expresses, among others, that there exist multiple sessions
%                               (say 1,2,3) for each subject -- i.e. the session property is *per
%                               subject*, and thus two sets tagged with session 3 do not
%                               necessarily refer to the same unique session (unless the subject,
%                               day, montage (if present) are also the same. If any such nesting
%                               relationship is properly expressed, the 'test_on' and 'per' settings
%                               will work as expected.
%
% Out:
%   Result : result of the partitioning; depends on the IndexSet
%
% See also:
%   set_partition, bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-29
dp;

args = arg_define(0:3, varargin, ...
    arg({'collection','Collection'},mandatory,[],'A dataset collection. Cell array of structs with meta-data.','type','expression'), ...
    arg({'inds','IndexSet'},mandatory,[],'Index set to use for partitioning. This partitioner is somewhat special: if this is [], the function generates the partitioning, i.e. returns a cell array of partitions, each of which is a cell array of {training collection, test collection, {''goal_identifier'',identifier for test collection}}. Otherwise, we assume that this is one of the partitioned collections that we generated, and just return it unmodified.','type','expression'), ...
    arg_sub({'settings','Settings'},{}, {
        arg({'scope_order','ScopeOrdering'},{'group','subject','day','montage','session','recording','block'},[],'Dataset identifiers ordered from coarsest to finest. This both defines what identifiers are known, and their hierarchical ordering.'), ...    
        arg({'restrict_to','RestrictTo'},[],[],'Restrict to data with these properties. Optionally a struct where one can set fields to values that all retained sets must have. E.g., by setting .day to 5, only sets with day=5 will be retained.','type','expression'), ...
        arg({'exclude','Exclude'},{},[],'Exclude data with these properties. Optionally a struct (or cell array of structs) where one can set fields to values that no retained set may have. Opposite of RestrictTo.'), ...
        arg({'test_on','TestOn'},'',[],'Separate test sets based on this property. This is a dataset identifier that determines which items to use as test sets. By default, this is the coarsest granularity (according to ScopeOrdering) that''s in the data.'), ...
        arg({'consider','Consider'},'each', {'each','all','last','allexceptfirst'}, 'Consider these items for testing. Allows to test on, for instance, each subject, or the last day of the experiment (using all others for training), or all days except the first (using the prior days for training).'), ...
        arg({'per','Per'},'',[], 'If consider is not each, determines the context to which the last/allexceptfirst refers (default: next coarsest granularity in the data according to ScopeOrdering).')...
    }, 'Settings for the partitioning.'));

[collection,inds,settings] = deal(args.collection, args.inds, args.settings);

% input validation
if ~iscell(collection) || ~all(cellfun('isclass',collection,'struct'))
    error('The given Collection argument must be a cell array of structs.'); end
if ~all(cellfun(@(s)isfield(s,'streams'),collection))
    disp_once('Note: each cell in the collection is expected to be a stream bundle (struct with .streams field), but not all are.'); end
if ~iscellstr(settings.scope_order)
    error('The given ScopeOrdering argument must be a cell array of strings, but was: %s',hlp_tostring(settings.scope_order)); end
if length(unique(settings.scope_order)) < length(settings.scope_order)
    error('The ScopeOrdering parameter must not contain duplicate elements, but was: %s',hlp_tostring(settings.scope_order)); end
if ~(isempty(settings.restrict_to) || isstruct(settings.restrict_to))
    error('If nonempty the given RestrictTo argument must be a struct, but was: %s',hlp_tostring(settings.restrict_to)); end
if ~(isempty(settings.exclude) || isstruct(settings.exclude))
    error('If nonempty the given Exclude argument must be a struct, but was: %s',hlp_tostring(settings.exclude)); end
if ~isempty(settings.test_on) && (~ischar(settings.test_on) || ~any(strcmp(settings.test_on,settings.scope_order)))
    error('The given TestOn setting should be empty or equal to one of the elements in ScopeOrdering (%s), but was: %s.',hlp_tostring(settings.scope_order),hlp_tostring(settings.test_on)); end
if ~any(strcmp(settings.consider,{'each','all','last','allexceptfirst'}))
    error('The given Consider option must be one of: {''each'',''all'',''last'',''allexceptfirst''}, but was: %s',hlp_tostring(settings.consider)); end
if ~isempty(settings.per) && (~ischar(settings.per) || ~any(strcmp(settings.per,settings.scope_order)))
    error('The given Per setting should be empty or equal to one of the elements in ScopeOrdering (%s), but was: %s.',hlp_tostring(settings.scope_order),hlp_tostring(settings.per)); end    

if ~isempty(inds)
    % check whether inds is of the right format
    if ~iscell(inds) || length(inds) ~= 3
        error('If the given indices are non-empty, the are expected to be a 3-element cell array, but were: %s',hlp_tostring(inds,1000)); end
    if ~all(cellfun('isclass',inds,'cell'))
        error('If the given indices are non-empty, it must be a cell array of cell arrays, but was: %s',hlp_tostring(inds,1000)); end
    % if so, pass it through
    res = inds;
    return;
end

known_granularities = settings.scope_order;

% first tag all items in the collection with a unique tracking id
% (will be used towards the end to determine the contents of the respective training sets)
for k=1:length(collection)
    collection{k}.tracking_index = k; end

% restrict the collection of test set material according to the restrict_to property
test_material = collection;
if ~isempty(settings.restrict_to)
    for fn = fieldnames(settings.restrict_to)'
        field = fn{1};
        present = find(cellfun(@(x)isfield(x,field),test_material));
        if isempty(present)
            fprintf('None of the sets in the collection has the property "%s"; ignoring this restriction.\n',field); 
        else
            test_material = test_material(present);
        end
        restrictor = settings.restrict_to.(field);
        retain = cellfun(isequal_weak(@(x)x.(field),restrictor),test_material,'UniformOutput',false);
        if ~any(retain)
            fprintf('None of the sets in the collection has value "%s" for property "%s"; ignoring this restriction.\n',hlp_tostring(restrictor),field); 
        else
            test_material = test_material(retain);
        end
    end
end

% implement the exclusion
if isstruct(settings.exclude)
    settings.exclude = {settings.exclude}; end
% for each excluder...
for e=1:length(settings.exclude)
    ex = settings.exclude{e};        
    % go through all fields and identify matches
    potential_match = true(1,length(test_material));
    for fn=fieldnames(ex)'
        field = fn{1};
        val = ex.(field);
        present = cellfun(@(x)isfield(x,field),test_material);
        matching_value = true(1,length(test_material));
        matching_value(present) = cellfun(@(x)isequal_weak(x.(field),val),test_material(present));
        potential_match = potential_match & present & matching_value;            
    end
    % retain only the non-matching test items
    test_material = test_material(~potential_match);
end

% remove all entries from known_granularities which are not contained in the test material and 
% remove all elements from the test material which are lacking one of the remaining granularity properties
retain = [];
for k=1:length(known_granularities)
    gran = known_granularities{k};
    present = cellfun(@(x)isfield(x,gran),test_material);
    if any(present)
        test_material = test_material(present);
        retain(end+1) = k; %#ok<*AGROW>
    end
end
known_granularities = known_granularities(retain);

if isempty(settings.test_on)
    % use coarsest granularity to test on
    if ~isempty(known_granularities)
        settings.test_on = known_granularities{1}; end
end

if ~strcmp(settings.consider,'each')         
    if isempty(settings.per)
        % need to determine the granularity at which to pack test items
        if isempty(settings.test_on)
            % we're testing on individual data items, so the finest known granularity in the 
            % data will be used for grouping
            if ~isempty(known_granularities)
                settings.per = known_granularities{end}; end
        else
            % see if we can find the next-coarser granularity
            known_gran = find(strcmp(settings.test_on,known_granularities),1);
            if ~isempty(known_gran)
                if known_gran == 1
                    error('The testing is already performed at the coarsest known granularity; cannot determine a coarser granularity to group data into. You may specify the ''per'' property manually.');
                else
                    settings.per = known_granularities{known_gran-1};
                end
            end
        end
        if isempty(settings.per)
            error('Please specify a granularity (i.e. a property) over which the ''consider'' clause should apply.'); end
    end
    % prune the test material appropriately
    present = find(cellfun(@(x)isfield(x,settings.per),test_material));
    if isempty(present)
        error('None of the sets in the collection has the ''per'' property "%s".\n',settings.per);
    else
        % exclude all items that are lacking the given field
        test_material = test_material(present);
    end
end


% build the test sets, which is a cell array of collections -- one per test set
test_sets = {};
if isempty(settings.test_on)
    % test at the granularity of data sets in the collection (finest)
    for k=1:length(test_material)
        test_sets{k} = test_material(k); end
else
    % test at the given granularity
    present = find(cellfun(@(x)isfield(x,settings.test_on),test_material));
    if isempty(present)
        error('None of the sets in the collection has the ''test_on'' property "%s".\n',settings.test_on);
    else
        % exclude all items that are lacking the given field
        test_material = test_material(present);
    end

    % we partition the test material into test sets by unique elements of the 'test_on' property
    % and all coarser properties, if any; next: build the list of these properties
    pos = find(strcmp(settings.test_on,known_granularities));
    if ~isempty(pos)
        partition_properties = known_granularities(1:pos);
    else
        partition_properties = {settings.test_on};
    end

    % now partition
    values = cellfun(@(x)getprops_as_string(x,partition_properties),test_material,'UniformOutput',false);
    uniquevals = unique(values);
    for v=1:length(uniquevals)
        matches = strcmp(values,uniquevals{v});
        test_sets{v} = test_material(matches);
    end
end

% optionally re-pack and/or prune the test sets according to the consider clause
if ~strcmp(settings.consider,'each')
    % we repack the test sets into groups by unique elements of the 'per' property
    % and all coarser properties, if any; next: build the list of these properties
    pos = find(strcmp(settings.per,known_granularities));
    if ~isempty(pos)
        partition_properties = known_granularities(1:pos);
    else
        partition_properties = {settings.per};
    end

    % sanity check: make sure that each test set has only one unique 'per' identifier
    for k=1:length(test_sets)
        if length(unique(cellfun(@(x)getprops_as_string(x,partition_properties),test_sets{k},'UniformOutput',false))) > 1
            error('You are trying to operate on your %s''s in groups per %s, but some of the %s''s consist of items in multiple groups. You can likely avoid this problem by specifying including both properties in the ''scope_order'' list at the appropriate place.'); end
    end

    % now partition by 'per'
    test_groups = {};
    values = cellfun(@(x)getprops_as_string(x{1},partition_properties),test_sets,'UniformOutput',false);
    uniquevals = unique(values);
    for v=1:length(uniquevals)
        matches = strcmp(values,uniquevals{v});
        test_groups{v} = test_sets(matches);
    end

    % and re-pack/prune to get the final test sets
    test_sets = {};
    switch settings.consider
        case 'all'
            % flatten the test groups
            for g=1:length(test_groups)
                test_sets{g} = [test_groups{g}{:}]; end
        case 'last'
            % take only the last set per group (= the one with the highest value for the test_on property)
            for g=1:length(test_groups)
                idx = argmax(cellfun(@(x)x{1}.(settings.test_on),test_groups{g}));
                test_sets{g} = [test_groups{g}{idx}];
            end
        case 'allexceptfirst'
            % take all except for the first set per group (= all except for those with the lowest value for the test_on property)
            for g=1:length(test_groups)
                idx = argmin(cellfun(@(x)x{1}.(settings.test_on),test_groups{g}));
                test_sets{g} = [test_groups{g}{setdiff(1:end,idx)}]; 
            end
        otherwise
            error('Unsupported ''consider'' clause specified.');
    end
end


% build the final partition
res = {};
for s=1:length(test_sets)
    testset = test_sets{s};
    % derive the corresponding training sets
    test_indices = cellfun(@(x)x.tracking_index,testset);
    not_in_test = cellfun(@(x)~ismember(x.tracking_index,test_indices),collection);
    trainset = collection(not_in_test);
    % get rid of tracking indices
    for k=1:length(trainset)
        trainset{k} = rmfield(trainset{k},'tracking_index'); end
    for k=1:length(testset)
        testset{k} = rmfield(testset{k},'tracking_index'); end
    % also derive the common identifying information of the test sets:
    common_fields = setdiff(fieldnames(testset{1}),{'streams'});
    for k=2:length(testset)
        common_fields = intersect(common_fields,fieldnames(testset{k})); end
    if length(testset) > 1
        % restrict to those fields that are unique across all test sets
        retain = [];
        for c=1:length(common_fields)
            field = common_fields{c};
            val = testset{1}.(field);
            if all(cellfun(@(x)isequal_weak(x.(field),val),testset(2:end)))
                retain(end+1) = c; end
        end
        common_fields = common_fields(retain);
    end
    % form an identifier struct for the test set; the cross-validation will append the contents
    % of this cell array as additional arguments into the trainer function
    identifier = [];
    for f=1:length(common_fields)
        identifier.(common_fields{f}) = testset{1}.(common_fields{f}); end
    % combine
    res{s} = {trainset,testset,{'goal_identifier',identifier}};
end


% obtain a string version of values from a struct, for a given cell array of properties
function y = getprops_as_string(x,props)
y = hlp_tostring(cellfun(@(p)x.(p),props,'UniformOutput',false));
