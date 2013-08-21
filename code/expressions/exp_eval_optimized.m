function varargout = exp_eval_optimized(x,iters)
% Evaluate the given expression using optimizations.
% Out-Args = exp_eval_optimized(exp)
%
% This function should produce the same result as exp_eval(), however it may cache intermediate
% results in memory or on disk (and reuse them later), or perform some (known-to-be-safe) reorderings
% of computations. Currently, this function is primarily relevant for evaluating filter expressions
% efficiently, and it contains some domain-specific optimizations to expedite their evaluation.
%
% In:
%   Expression : expression to evaluate
%
%   Iterations : optionally restrict the number of iterations done by exp_eval (default: 1)
%
% Out:
%   Out-Args   : the result of the evaluation
%
% Notes:
%   The currently implemented optimizations are memory and optional disk caching of all expressions
%   (while it does not attempt to disk-cache expressions that contain a set_partition() since, as 
%   there can be 100s of partitioned versions during bootstrap/resampling calculations). The other
%   optimization is that partitioning steps are being reodered to come after expensive filter operations
%   (e.g., resampling), so that the same data does not have to be processed over and over for different
%   overlapping partitions. Partitioning steps are not moved over processing steps that compute
%   statistics over the data (such as ICA, SSA or standardization), as declared by the 
%   'independent_trials' property of the respective stage.
%
% See also:
%   exp_eval
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15

if ~exist('iters','var')
    iters = 1; end

varargout = {x};
if isfield(x,{'head','parts'})
    % first (and key) optimization: "pull up" any set_partition() node in the expression as high as 
    % possible (i.e. partition the data as late as possible); the only nodes that cannot be skipped by 
    % a partition operation are those which have their 'independent_trials' property set to false.
    x = pull_up(x);

    % now add caching hints recursively...
    [x,unpartitioned] = add_cache_hints(x,iters);
    
    % ... and evaluate with the remaining evaluation hints around the top-level expression
    [varargout{1:nargout}] = hlp_scope({'memoize',{'memory',1,'disk',double(unpartitioned)}}, @exp_eval,x,iters);
end


% add disk cache hints around whatever is sitting below a set_partition expression...
function [x,unpartitioned] = add_cache_hints(x,iters)
unpartitioned = true;
if all(isfield(x,{'parts','head'}))
    % first recurse and find out whether the sub-expressions are partitioned or not
    for p=1:length(x.parts)
        [x.parts{p},unp] = add_cache_hints(x.parts{p},iters);
        unpartitioned = unpartitioned & unp;
    end

    % check if this is a set_partition expression
    chead = char(x.head);
    if strcmp(chead,'set_partition')
        if unpartitioned
            % it is sitting on top of an unpartitioned sub-expression: add a disk cache hint
            % around that sub-expression
            x.parts{1} = exp_block({exp_rule(@memoize,{'memory',1,'disk',1})},x.parts{1},iters);
        end
        unpartitioned = false;
    end
end


function [x,chead] = pull_up(x)
% recursively pull up set_partition(x,_) nodes in an expression structure
chead = char(x.head);
num_expressions = 0;                                    % number of sub-expressions
partition_at = [];                                      % indices of the set_partition sub-expressions in parts, if any

% recurse into parts...
for p=1:length(x.parts)
    if isfield(x.parts{p},{'head','parts'})
        [x.parts{p},head_p] = pull_up(x.parts{p});
        % count the number of sub-expressions
        num_expressions = num_expressions+1;
        % remember the places where we have set_partition sub-expressions
        if strcmp(head_p,'set_partition')
            partition_at(end+1) = p; end
    end
end

% now perform the reordering, if applicable
if num_expressions == length(partition_at) && num_expressions > 0
    % check if the current node has the independent_trials=false property (which prevents reordering)
    props = arg_report('properties',x.head);
    if ~props.independent_trials
        return;  end
    
    % we can proceed
    if num_expressions == 1
        % simple standard case: reorder this node with the set_partition below it
        % make a backup of our old set_partition(<y>,index_set) sub-expression
        oldchild = x.parts{partition_at};
        % replace our old set_partition(<y>,index_set) by the <y>
        x.parts{partition_at} = oldchild.parts{1};
        % and replace the whole thing by a new set_partition(<the whole thing>,index_set) node
        x = struct('head',{@set_partition},'parts',{{x,oldchild.parts{2:end}}});
        chead = 'set_partition';
        
        % finally recurse again and pull up any subordinate set_partition's
        for p=1:length(x.parts)
            if isfield(x.parts{p},{'head','parts'})
                x.parts{p} = pull_up(x.parts{p}); end
        end
    else
        % fusion of multiple partitioned data sets (only permitted if indices are identical)
        return; % not yet implemented...
    end
end
