function res = is_atom(x)
% Determine whether an expression can be further subdivided into parts
% literal matlab expressions (including matrices, most structs, strings, function_handles) are considered atomic
%
% The exceptions are:
%  * struct with .head / .part fields      (canonically represented expression)
%  * struct with .tracking.expression field  (impure expression, i.e. expression with value)
%
% Notes:
%  in the future, more data structures may be made accessible to the expression system
%
% See also:
%   exp_eval, exp_head, exp_parts
%
%                                                     Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                                     2010-04-15

if isfield(x,{'head','parts'})
    res = true;
elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
    res = is_atom(x.tracking.expression);
else
    res = false;
end