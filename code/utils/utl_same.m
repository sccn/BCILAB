function r = utl_same(a,b)
% Check whether two expressions are structurally identical.
%
% Differs from isequal() in the following ways:
%  a) ignores the values associated with impure expressions
%  b) compares lambda expressions by their string representation; necessary for exp_symbol to work as expected
%  note: exp_symbol('x') and @x are considered non-equal by utl_same.
%
% Note: This function treats anonymous functions nested inside cell arrays or structs as non-equal,
%       unless they are copies of each other; use isequal_weak for this.
%
% See also:
%   exp_same
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-04-19

if isstruct(a) && isstruct(b)
    % structures may be impure expressions    
    if isfield(a,'tracking') && isfield(a.tracking,'expression')
        a = a.tracking.expression; end
    if isfield(b,'tracking') && isfield(b.tracking,'expression')
        b = b.tracking.expression; end
    
    % or canonically represented expressions
    if isfield(a,{'head','parts'}) & isfield(b,{'head','parts'}) %#ok<AND2>
        if ~isequal(a.head,b.head) || (length(a.parts) ~= length(b.parts))
            r = false;
        else
            for k=1:length(a.parts)
                if ~utl_same(a.parts{k},b.parts{k})
                    r = false;
                    return;
                end
            end
            r = true;
        end
    else
        % or general structures
        r = isequalwithequalnans(a,b);
    end        
elseif isa(a,'function_handle') && isa(b,'function_handle')
	% function handles are compared via their string representation
    r = strcmp(char(a),char(b));
else
    % everything else is compared using isqual
    r = isequalwithequalnans(a,b);
end
