function r = is_evaluated(x)
% determine whether an expression can be taken as a fully evaluated MATLAB value
% this holds for both impure expressions and literal expressions (almost all MATLAB data structures; see exp_eval for a detailed list)
% it excludes:
%   * canonical expression structs (with the only fields being head and parts)
%   * undefined functions (which are not valid MATLAB values and serve as a shorthand notation for unevaluated symbols, e.g. @x)
%   * lambda functions which yield canonical expression structs when invoked (such as those lambdas created by exp_symbol or exp_declare)
%
% Notes:
%   Evaluated(x) := impure(x) || ~(canonical(x) || undefined(x) || symbolic_lambda(x));
%
% See also:
%   exp_eval

r = (isfield(x,'tracking') && isfield(x.tracking,'expression')) || ~( all(isfield(x,{'head','parts'})) || (isa(x,'function_handle') && (is_undefined_function(x) || ~isempty(regexp(char(x),'@\(varargin\)struct\(''head'',\{.*\},''parts'',\{varargin\}\)','once')))));