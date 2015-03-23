function targets = utl_collection_targets(testcollection)
% Internal. Obtain the target values for a dataset collection, as part of a cross-validation.
% Targets = utl_collection_targets(TestCollection)
%
% In:
%   TestCollection : dataset collection (cell array) to which a model shall be applied; the elements
%                    can be stream bundles or EEGLAB dataset structs
%
% Out:
%   Targets : extracted target values for each collection element, concatenated
%
% See also:
%   utl_collection_tester
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-29
dp;

% input validation
if ~iscell(testcollection) || ~all(cellfun('isclass',testcollection,'struct'))
    error('The given TestCollection argument must be a cell array of structs, but was: %s',hlp_tostring(testcollection,10000)); end

% note: these are drawn from the first stream (as all streams have the same markers)
for k=length(testcollection):-1:1
    if isfield(testcollection{k},'streams') && iscell(testcollection{k}.streams) && ~isempty(testcollection{k}.streams)
        targets{k} = set_gettarget(testcollection{k}.streams{1}); 
    else
        targets{k} = set_gettarget(testcollection{k}); 
    end
end
targets = utl_aggregate_results(targets{:});
