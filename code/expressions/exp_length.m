function l = exp_length(x)
% Get the number of elements (parameters) in any expression.
% Length = exp_length(Expression)
%
% In:
%   Expression  : any expression
%
% Out:
%   Length      : the expression's length; highest allowed index allowed for exp_part.
%
% See also:
%   exp_head, exp_part, exp_parts
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-26

if isfield(x,{'head','parts'})
    % struct with fields .head & .parts
    l = length(x.parts);
elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
    % impure expressions viewed as expressions
    l = exp_length(x.tracking.expression);
else
    % most MATLAB data structures are not further broken down when reflected into the expression 
    % system
    l = 1;
end