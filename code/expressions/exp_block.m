function res = exp_block(rules, expr, iters) 
% Define a block with dynamically scoped values for symbols used in the contained expression.
% Result = exp_block(Rules, Expression, Iterations)
%
% In:
%   Rules : cell array of exp_set/exp_setdelayed/exp_rule expressions; the left-hand side
%           (first argument) of the assignment expressions is replaced by the respective
%           right-hand side (second argument) in the scope of the block expression; this
%           implements dynamically scoped variables.
%
%   Expression  : an expression to be evaluated with assignments as listed in the block
%
%   Iterations : optionally set the maximum number of iterations done by exp_eval (default: 1)
%                for real symbolic computations, use Inf here
%
% Out:
%   Result : a Block expression
%
% See also:
%   hlp_scope, exp_rule
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-03


if ~exp_beginfun('symbolic') return; end

if ~exist('iters','var')
    iters = 1; end
if ~iscell(rules)
    error('Rules should be given as a cell array of the form {Set(...),Set(...),SetDelayed(...),Set(...), ...}'); end

% turn the assignments into a struct
assignments = struct();
for r=1:length(rules)
    rule = rules{r};
    if any(strcmp(char(rule.head),{'exp_set','exp_rule'}))
        rhs = exp_eval(rule.parts{2},inf);
    else
        rhs = rule.parts{2};
    end    
    assignments.(char(get_function_symbol(rule.parts{1}))) = rhs;    
end

% now evaluate the expression with assignments applied using hlp_scope
res = hlp_scope(assignments,@exp_eval,expr,iters);

exp_endfun;
