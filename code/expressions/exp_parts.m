function p = exp_parts(x)
% Get the parts (parameters) of any expression.
% Parts = exp_parts(Expression)
% 
% In:
%   Expression  : any expression
%
% Out:
%   Parts       : The parts of the expression, a cell array of expressions. An expression written as
%                 f(x,y,z) has {x,y,z} as its parts (though the data structure by which the
%                 expression is represented in MATLAB is a (canonical expression) struct). Since
%                 'sdfsf' and 10.3 are also expressions, their head can be obtained, too.
%                  * the parts of a canonical expression struct s are s.parts
%                  * the parts of an impure expression struct s are the parts of
%                    s.tracking.expression.
%                  * the parts of a native MATLAB data structure depend on the data type:
%                    - single-row cell array: the cell array itself 
%                    - any other value x: {x) (this is a catch-all, which may later be further 
%                      refined)
%
% See also:
%   exp_head, exp_part
%
% Note:
%   There is one type of function_handle that receives special treatment, namely the lambda
%   functions created by exp_symbol; their parts is the contained MATLAB function handle (what was
%   passed to exp_symbol, with an @ prepended), wrapped in a cell array.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

while isfield(x,'tracking') && isfield(x.tracking,'expression')
    x = x.tracking.expression; end

if isfield(x,{'head','parts'})
    % canonically represented expression, struct with fields .head * .parts
    p = x.parts;
elseif has_function_symbol(x)
    % function handle which is reflected as a symbol, e.g. @name
    p = {get_function_symbol(x)};
else
    % anything else
    p = {x};
end