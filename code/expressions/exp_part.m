function p = exp_part(x,idx)
% Get a specific part of an expression; can be more efficient than exp_parts.
% Parts = exp_part(Expression, Index)
% 
% In:
%   Expression  : any expression
%   Index       : index of the part; 
%                  * 0 is the expression's head
%                  * 1..n is the expression's parts, n being exp_length(Expression)
%                  * negative indices count parts from the end
%                  * index vectors denote paths into the expression tree
%
% Out:
%   Part        : the selected part of the expression, as in exp_parts
%
% See also:
%   exp_head, exp_parts, exp_length
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-26

while isfield(x,'tracking') && isfield(x.tracking,'expression')
    x = x.tracking.expression; end

try
    if length(idx)==1
        if idx > 0
            % positive indices behave as exp_parts(x){idx} (if that could be written in MATLAB...)
            if isfield(x,{'head','parts'})
                p = x.parts{idx};
            elseif has_function_symbol(x)
                p = get_function_symbol(x);
            elseif idx == 1
                p = x;
            else
                error; %#ok<LTARG>
            end
        elseif idx == 0
            % part 0 is the head
            p = exp_head(x);
        elseif idx < 0
            % negative indices count from the end (-1 is the last part)
            p = exp_part(x,exp_length(x)+idx+1);
        end
    else
        % index vectors specify paths into the expression
        p = exp_part(exp_part(x,idx(1)),idx(2:end));
    end
catch
    error('out-of-range part accessed.');
end