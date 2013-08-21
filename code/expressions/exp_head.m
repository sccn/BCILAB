function h = exp_head(x)
% Get the head (=function label) of any expression.
% Head = exp_head(Expression)
% 
% In:
%   Expression  : any expression
%
% Out:
%   Head : The head of the expression, again an expression. An expression written as 
%          f(x,y,z) has f as its head (though the data structure by which the expression is
%          represented in MATLAB is a (canonical expression) struct). Since {x,y,z} and 10.3 are
%          also expressions, their head can be obtained, too.
%           * the head of a canonical expression struct s is s.head
%           * the head of an impure expression struct s is the head of s.tracking.expression.
%           * the head of a native MATLAB data structure depends on the data type:
%             - scalar number: @Number
%             - single-row char array: @String
%             - function_handle: @Symbol
%             - anything else: @Value (this is a catch-all, which may later be further refined into 
%                              things like @Struct, etc.)
%
% See also:
%   exp_part, exp_parts, exp_eval, exp_eval_optimized
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

while isfield(x,'tracking') && isfield(x.tracking,'expression')
    x = x.tracking.expression; end

if isfield(x,'head')
    % structure with fields .head & .parts, as created by most of the expression machinery
    h = x.head;
elseif isa(x,'function_handle')
    % all functions are considered symbols, even though some of them may be lambda functions
    h = @Symbol;
elseif ischar(x) && size(x,1) == 1
    % regular strings; Nd char arrays are currently reflected as Value (rather than
    % List(String(...), String(...), ...))
    h = @String;
elseif isnumeric(x) && isscalar(x)
    % regular numbers; matrices are currently reflected as Value (rather than 
    % List(List(Number, Number, ...), ...))
    h = @Number;
else
    % generic value
    h = @Value;
end