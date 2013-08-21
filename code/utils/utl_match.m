function [r,dict] = utl_match(x,f,dict)
% Check whether an expression is matched by a pattern. 
% [Equals,Assignment] = utl_match(Expression, Pattern, Assignment)
%
% The pattern may contain blanks, conditional expressions, and named sub-patterns.
% Assignment is an optional output which contains sub-expressions matched by named parts of the pattern.
%  note: exp_symbol('x') and @x are considered non-equal by utl_match.
%
% See also:
%   exp_match, utl_same
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-04-19

if isfield(f,{'head','parts'})
    % f is a canonical expression (may be a pattern or a blank)
    head = char(f.head);
    if strcmp(head,'Pattern')
        % f is a pattern: try to match its form
        [r,dict] = utl_match(x,f.parts{2},dict);
        if r
            % we have a match for a named pattern
            name = char(f.parts{1});
            if name(1) == '@'
                % normalize if it has a symbolic lambda name
                name = name(28:end-21); end
            % try to store the name-x relationship in the dictionary
            [r,dict] = update_dict(dict,name,x);
        end
    elseif strcmp(head,'Blank') && (isempty(f.parts) || strcmp(char(get_function_symbol(f.parts{1})),char(exp_head(x))))
        % we have a match for an unnamed pattern
        r = 1;
    elseif strcmp(head,'Condition')
        % f is a condition expression: try to match its expression
        [r,dict] = utl_match(x,f.parts{1},dict);
        % then substitute the named sub-expressions into f's condition and evaluate that
        r = r && exp_eval(utl_replacerepeated(f.parts{2},dict),inf);
    else
        % f is a non-pattern expression... match recursively
        if isfield(x,'tracking') && isfield(x.tracking,'expression')
            % normalize impure x
            x = x.tracking.expression; end
        if isfield(x,{'head','parts'})
            % x has substructure: we have to match recursively
            if ~isequal(x.head,f.head) || (length(x.parts) ~= length(f.parts))
                r = false;
            else
                for k=1:length(x.parts)
                    [r,dict] = utl_match(x.parts{k},f.parts{k},dict);
                    if ~r
                        return; end
                end
                r = true;
            end
        else
            % x has no substructure (but f has)
            r = false;
        end        
    end
elseif isfield(f,'tracking') && isfield(f.tracking,'expression')
    % f is an impure expression; descend
     [r,dict] = utl_match(x,f.tracking.expression,dict);
elseif isa(f,'function_handle') && isa(x,'function_handle')
    % f an x are function handles: compare their string representations
    r = strcmp(char(f),char(x));
else
    % f is atomic; x must be identical
    r = isequalwithequalnans(f,x);
end

% update the given dictionary with a name-value pair; also return whether the naming is still consistent
function [ok,dict] = update_dict(dict,name,val)
ok = true;
if isfield(dict,name)
    if ~utl_same(dict.(name),val)
        % mismatch!
        ok = false;
        return;
    end
else
    dict.(name) = val;
end
