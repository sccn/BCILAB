function res = is_symbolic_lambda(x)
% internal: a symbolic lambda function is one which generates expressions when invoked with arguments (this is what exp_symbol generates)
res = isa(x,'function_handle') && ~isempty(regexp(char(x),'@\(varargin\)struct\(''head'',\{.*\},''parts'',\{varargin\}\)','once'));