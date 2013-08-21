function x = utl_purify_expression(x)
% Recursively remove impure elements from an expression. 
% Output-Expression = utl_purify_expression(Input-Expression)
%
% In particular, x.tracking.expression is replaced by x. It is assumed that the head is not an impure expression.
%
%                                                     Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                                     2010-05-23

if isfield(x,{'head','parts'})
    % x is a canonical expression: recurse
    for i=1:length(x.parts)
        x.parts{i} = utl_purify_expression(x.parts{i}); end
elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
    % x is an impure expression: purify & recurse
    x = utl_purify_expression(x.tracking.expression);
end
