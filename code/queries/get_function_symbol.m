function result___ = get_function_symbol(expression___)
% internal: some function_handle expressions have a function symbol (an @name expression), and this function obtains it

% note: we are using funny names here to bypass potential name conflicts within the eval() clause further below
if ~isa(expression___,'function_handle')
    error('the expression has no associated function symbol.'); end

string___ = char(expression___);
if string___(1) == '@'
    % we are dealing with a lambda function
    if is_symbolic_lambda(expression___)
        result___ = eval(string___(27:end-21));        
    else
        error('cannot derive a function symbol from a non-symbolic lambda function.'); 
    end
else
    % we are dealing with a regular function handle
    result___ = expression___;
end