function x = utl_releasehold(x)
% Peel off a layer of hold expressions from some expression.
% [Out-Expression] = utl_releasehold(In-Expression)
%
% In:
%   In-Expression      : some expression
%
% Out:
%   Out-Expression     : In-Expression with the first level of Hold's peeled off
%
% Notes:
%   Descends into cell arrays.
%
% See also:
%   exp_releasehold, exp_hold
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-05-19

if isfield(x,{'head','parts'})
    if strcmp(char(x.head),'Hold')
        % peel off a layer of hold (but do not recurse!)
        x = x.parts{1};
    else    
        % x is a canonical expression: recurse
        x.head = utl_releasehold(x.head);
        for i=1:length(x.parts)
            x.parts{i} = utl_releasehold(x.parts{i}); end
    end
elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
    % x is an impure expression: descend
    x.tracking.expression = utl_releasehold(x.tracking.expression);
elseif iscell(x) && size(x,1) == 1
    % x is a cell array
    for i=1:length(x)
        x{i} = utl_releasehold(x{i}); end
end
