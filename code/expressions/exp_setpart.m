function x = exp_setpart(x,idx,v)
% Set a specific part of an expression.
% Result = exp_setpart(Expression, Index, Value)
% 
% In:
%   Expression  : any expression
%   Index       : index of the part; 
%                  * 0 is the expression's head
%                  * 1..n is the expression's parts, n being exp_length(Expression)
%                  * negative indices count parts from the end
%                  * index vectors denote paths into the expression tree
%   Value       : the value to be assigned to the expression's part
%
% Out:
%   Result      : the expression with the referenced part replaced by Value
%
% See also:
%   exp_part, exp_length
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-09

try
    if length(idx)==1
        if idx > 0
            % regular indexing, analogous to exp_part
            if isfield(x,{'head','parts'})
                x.parts{idx} = v;
            elseif isfield(x,'tracking') && isfield(x.tracking,'expression')
                x.tracking.expression = exp_setpart(x.tracking.expression,idx,v);
            elseif has_function_symbol(x) 
                if ~isequal(exp_part(x,idx),v)
                    error('Cannot assign to a symbol.'); end
            elseif idx == 1
                x = v;
            else
                error; %#ok<LTARG>
            end
        elseif idx == 0
            % assignment of the head
            if isfield(x,'head')
                if ~isequal(exp_head(x),v) && any(cellfun(@(b)utl_same(v,b),{@Symbol,@Number,@Value,@List,@String}))
                    % turning a canonical expression into a builtin one
                    error('Head assignments that change the representation of an expression are currently not allowed.'); end
                x.head = v;
            elseif ~isequal(exp_head(x),v)
                % turning a builtin expression into a canonically represented one 
                error('Head assignments that change the representation of an expression are currently not allowed.'); 
            end
        elseif idx < 0
            % negative indices count from the end (-1 is the last part)            
            x = exp_setpart(x,exp_length(x)+idx+1,v);
        end
    else
        % index vectors specify paths into the expression
        x = exp_setpart(x,idx(1),exp_setpart(exp_part(x,idx(1)),idx(2:end),v));
    end
catch
    error('invalid access to a part.');
end