function res = utl_collection_partition(varargin)
% Partition a dataset collection according to some settings
% Result = utl_collection_partition(Collection,IndexSet,Settings)
%
% This function is a valid partitioner for utl_crossval, utl_nested_crossval, and utl_searchmodel,
% which accepts a collection (cell array) of datasets. The datasets may have additional meta-data
% fields that describe, for instance, the subject id ('subject'), session number ('session'), etc.,
% of the data. For the meta-data field that are recognized by default, see default value of the
% argument ScopeOrdering (which can also be overridden).
%
% The settings of this function (which are accessible in bci_train as the EvaluationScheme)
% basically allow one to set up how such a dataset collection shall be partitioned in a
% cross-validation.
%
% In:
%   Collection : A dataset collection, i.e. cell array of structs (each of which may have meta-data)
%
%   IndexSet : The index set to use for partitioning; this partitioner is somewhat special;
%              * If this is [], the function generates the partitioning, i.e. returns
%                a cell array of partitions, each of which is a cell array of:
%                {training-collection, test-collection, {'goal_identifier',identifier for test collection}}
%                (the thrid cell array is used by utl_crossval as control arguments for the training
%                process that receives the training collection)
%
%              * If this is not [], we assume that this is one of the partitioned collections that 
%                we generated (e.g., training-collection), and just return it unmodified.
%
%   Settings : Settings that govern how a collection of recordings shall be partitioned into folds 
%              and train/test sets for cross-validation:
%
%               RestrictTestsetsTo : A struct that may have fields and corresponding values that
%                                    allow to restrict what recordings shall be allowed in test sets
%                                    (only recordings that have each of the given fields and
%                                    corresponding values will be used for testing). All recordings
%                                    not used for testing are used for training. For instance, one
%                                    may choose to generally only test on recordings with 'session'
%                                    set to 'driving'. (default: blank / no restriction)
%               ExcludeFromTestsets : A struct that may have a combination of fields and corresponding
%                                    values that allow to restrict what recordings shall be used for
%                                    testing (recordings that have all of the given fields and
%                                    corresponding values will be excluded from testing). It is also
%                                    possible to give a cell array of structs to exclude multiple
%                                    combinations. All recordings not used for testing are used for
%                                    training. For instance, one may choose to exclude certain kinds
%                                    of training sessions from testing. (default: blank / no
%                                    restriction)
%               TestUnit : Group test recordings into units according to this property. The retained
%                          test recordings are grouped together into units such that each unit has a
%                          unique value for this property (e.g., if this is set to ''day'', then all
%                          datasets for the given day are a considered a single test unit). In the
%                          simplest case, one test unit corresponds to one test set during
%                          cross-validation. As a result, the cross-validation would have as many
%                          folds as there were days in the dataset collection. By default this is
%                          the coarsest property that''s present in the data (e.g., subject).
%               Consider : A modifier that determines how to pack and optionally prune the test 
%                          units into final test set used in the cross-validation. If this is
%                          'each', then each test unit corresponds to one test fold. Otherwise, this
%                          argument and the Per argument (which is a property name like TestUnit)
%                          forms a description of how to group the test sets used in the
%                          cross-validation, and reads (as a template): "Consider <<Consider>>
%                          <<TestUnit>> per <<Per>> as a test set." -- e.g., "Consider 'last' 'day'
%                          per 'subject' as test set." when Consider='last', TestUnit='day', and
%                          Per='subject'. This means that for each subject, the recordings taken on
%                          the last day go into a single test set (and there is one such test set per
%                          subject). The allowed values are: 'each', 'all' (all test units per X
%                          form a single test set), 'last' (as explained), 'allexceptfirst' (all
%                          except the first test units per X form a single test set). (default: 'each')
%               Per : This is ignored if Consider is 'each' (the default). Otherwise, this is the
%                     property according to which test data are separated into folds (so, if this is 
%                     'session' then there will be one test set per session in the data) What goes
%                     into that fold depends on Consider and TestUnit (e.g., the last block of the
%                     session, or all except the first blocks, if Consider is 'allexceptfirst' and
%                     TestUnit is 'block'). (default: next coarser granularity in the data after 
%                     TestUnit)
%               ScopeOrdering : Odering of partitioning-relevant properties from coarsest granularity 
%                               to finest (default: {'group', 'subject', 'day', 'montage',
%                               'session', 'recording','block'}) This expresses, among others, that
%                               there exist multiple sessions (say 1,2,3) for each subject -- i.e.
%                               the session property is *per subject*, and thus two sets tagged with
%                               session 3 do not necessarily refer to the same unique session
%                               (unless the subject, day, montage (if present) are also the same. If
%                               any such nesting relationship is properly expressed, the TestUnit
%                               and Per settings will work as expected.
%               NoCrossvalidation : Whether to disable the cross-validation (like passing 0 to
%                                   bci_train)
%
% Out:
%   Result : result of the partitioning; depends on the IndexSet
%            * If the IndexSet is [], then this is a cell array of cell arrays, each of which has
%              the form {train-partition, test-partition, trainer-arguments} where *-partition is a
%              collection (cell array) of dataset structs, and trainer-arguments is a cell array of
%              arguments that are meant to be passed into the training function that processes the
%              training collection (this is done by utl_crossval).
%            * If the IndexSet is a cell array of dataset structs, it will be passed right through
%              unmodified.
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
        arg({'restrict_to','RestrictTestsetsTo'},[],[],'Restrict test sets to these properties. Optionally a struct where one can set fields to values that all retained sets must have. E.g., by setting .day to 5, only sets with day=5 will be retained.','type','expression'), ...
        arg({'exclude','ExcludeFromTestsets'},{},[],'Exclude test sets with these properties. Optionally a struct (or cell array of structs) where one can set fields to values that no retained set may have. Opposite of RestrictTo.','type','expression'), ...
        arg({'test_unit','TestUnit'},'',[],'Group test sets into units according to this property. The remaining test sets are grouped together into units such that each unit has the same value of this property (e.g., if this is set to ''day'', then all datasets for the given day are a considered a single test set. As a result, the cross-validation will have as many folds as there are day in the collection. By default this is the coarsest property that''s in the data (e.g., subject).'), ...
        arg({'consider','Consider'},'each', {'each','all','last','allexceptfirst'}, 'Consider these items for testing. Allows to test on, for instance, each subject, or the last day of the experiment (using all others for training), or all days except the first (using the prior days for training).'), ...
        arg({'per','Per'},'',[], 'Split into one fold per. If consider is not each, determines the context to which the last/allexceptfirst refers (default: next coarsest granularity in the data according to ScopeOrdering).')...
        arg({'scope_order','ScopeOrdering'},{'group','subject','day','montage','session','recording','block'},[],'Dataset identifiers ordered from coarsest to finest. This both defines what identifiers are known, and their hierarchical ordering.'), ...    
        arg({'no_crossval','NoCrossvalidation'},false,[],'Whether to disable the cross-validation. Evaluation measures will be NaN.'), ...    
    }, 'Settings for the partitioning.'));

[collection,inds,settings] = deal(args.collection, args.inds, args.settings);

if settings.no_crossval && isempty(inds)
    res = {{collection, {}, {}}}; 
    return;
end

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
if ~isempty(settings.test_unit) && (~ischar(settings.test_unit) || ~any(strcmp(settings.test_unit,settings.scope_order)))
    error('The given TestUnit setting should be empty or equal to one of the elements in ScopeOrdering (%s), but was: %s.',hlp_tostring(settings.scope_order),hlp_tostring(settings.test_unit)); end
if ~any(strcmp(settings.consider,{'each','all','last','allexceptfirst'}))
    error('The given Consider option must be one of: {''each'',''all'',''last'',''allexceptfirst''}, but was: %s',hlp_tostring(settings.consider)); end
if ~isempty(settings.per) && (~ischar(settings.per) || ~any(strcmp(settings.per,settings.scope_order)))
    error('The given Per setting should be empty or equal to one of the elements in ScopeOrdering (%s), but was: %s.',hlp_tostring(settings.scope_order),hlp_tostring(settings.per)); end    

if ~isempty(inds)
    % check whether inds is of the right format
    if ~iscell(inds) 
        error('If the given indices are non-empty, the are expected to be a cell array, but were: %s',hlp_tostring(inds,1000)); end
    if ~all(cellfun('isclass',inds,'struct'))
        error('If the given indices are non-empty, it must be a cell array of structs, but was: %s',hlp_tostring(inds,1000)); end
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
    if ~isempty(fieldnames(ex))
        % mask of items that might be matched by this excluder
        potential_match = true(1,length(test_material));
        % go through all fields and identify matches
        for fn=fieldnames(ex)'
            field = fn{1};
            val = ex.(field);
            present = cellfun(@(x)isfield(x,field),test_material);
            matching_value = true(1,length(test_material));
            matching_value(present) = cellfun(@(x)isequal_weak(x.(field),val),test_material(present));
            % remove datasets that don't have the respective excluder field or value from the matching set
            potential_match = potential_match & present & matching_value;            
        end
        % retain only the non-matching test items
        test_material = test_material(~potential_match);
    end
end

% remove all entries from known_granularities which are not contained in the candidate test material
% and remove all elements from the candidate test material which are lacking one of the remaining
% granularity properties that other test candidates have
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

if isempty(settings.test_unit)
    % use coarsest granularity that's in the data to test on
    if ~isempty(known_granularities)
        settings.test_unit = known_granularities{1}; end
end

% in anything but the trivial 'each' case we need to resolve the 'per' property fully, 
% (and should also remove any test sets that don't have that property)
if ~strcmp(settings.consider,'each')         
    if isempty(settings.per)
        % need to determine the granularity at which to pack test items
        if isempty(settings.test_unit)
            % we're testing on individual data items, so the finest known granularity in the 
            % data will be used for grouping
            if ~isempty(known_granularities)
                settings.per = known_granularities{end}; end
        else
            % see if we can find the next-coarser granularity
            known_gran = find(strcmp(settings.test_unit,known_granularities),1);
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
if isempty(settings.test_unit)
    % all we can do is test at the granularity of data sets in the collection (finest)
    for k=1:length(test_material)
        test_sets{k} = test_material(k); end
else
    % we group the set sets by unique values of the "test_unit" property
    present = find(cellfun(@(x)isfield(x,settings.test_unit),test_material));
    if isempty(present)
        error('None of the sets in the collection has the ''test_unit'' property "%s".\n',settings.test_unit);
    else
        % exclude all items that are lacking the given field
        test_material = test_material(present);
    end

    % we partition the test material into test sets by unique elements of the 'test_unit' property
    % and all coarser properties, if any; next: build the list of these properties
    pos = find(strcmp(settings.test_unit,known_granularities));
    if ~isempty(pos)
        partition_properties = known_granularities(1:pos);
    else
        partition_properties = {settings.test_unit};
    end

    % now partition
    values = cellfun(@(x)getprops_as_string(x,partition_properties),test_material,'UniformOutput',false);
    uniquevals = unique(values);
    for v=1:length(uniquevals)
        matches = strcmp(values,uniquevals{v});
        test_sets{v} = test_material(matches);
    end
end

% finally (optionally) re-group and/or prune the test sets according to the consider clause
if ~strcmp(settings.consider,'each')
    % we regroup the test sets by unique values of the 'per' property and all coarser properties,
    % if any
    
    % build the list of these properties
    pos = find(strcmp(settings.per,known_granularities));
    if ~isempty(pos)
        partition_properties = known_granularities(1:pos);
    else
        partition_properties = {settings.per};
    end

    % sanity check: make sure that each test set has only one unique 'per' identifier
    % otherwise we'd have to break up the units (if you get this error then it's possible that your
    % ScopeOrdering argument does not correctly respect the natural scopes in your data)
    for k=1:length(test_sets)
        if length(unique(cellfun(@(x)getprops_as_string(x,partition_properties),test_sets{k},'UniformOutput',false))) > 1
            error('You are trying to operate on your %s''s in groups per %s, but some of the %s''s consist of items in multiple groups. You can likely avoid this problem by specifying including both properties in the ''scope_order'' list at the appropriate place.'); end
    end

    % now partition the test sets by unique values of 'per'
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
            % take only the last set per group (= the one with the highest value for the test_unit property)
            for g=1:length(test_groups)
                idx = argmax(cellfun(@(x)x{1}.(settings.test_unit),test_groups{g}));
                test_sets{g} = [test_groups{g}{idx}];
            end
        case 'allexceptfirst'
            % take all except for the first set per group (= all except for those with the lowest value for the test_unit property)
            for g=1:length(test_groups)
                idx = argmin(cellfun(@(x)x{1}.(settings.test_unit),test_groups{g}));
                test_sets{g} = [test_groups{g}{setdiff(1:end,idx)}]; 
            end
        otherwise
            error('Unsupported ''consider'' clause specified.');
    end
end

% get rid of tracking indices again
collection_notracking = collection;
for k=1:length(collection_notracking)
    collection_notracking{k} = rmfield(collection_notracking{k},{'tracking_index'}); end

% build the final partition
res = {};
for s=1:length(test_sets)
    testset = test_sets{s};
    % derive the corresponding training sets
    test_indices = cellfun(@(x)x.tracking_index,testset);
    not_in_test = cellfun(@(x)~ismember(x.tracking_index,test_indices),collection);
    trainset = collection_notracking(not_in_test);
    % get rid of tracking indices
    for k=1:length(testset)
        testset{k} = rmfield(testset{k},{'tracking_index'}); end
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

1; % debug breakpoint


% obtain a string version of values from a struct, for a given cell array of properties
function y = getprops_as_string(x,props)
y = hlp_tostring(cellfun(@(p)x.(p),props,'UniformOutput',false));
