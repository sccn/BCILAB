function res = is_builtin_function(x)
% check if the given object is a valid handle to a builtin function
res = isa(x,'function_handle') && exist(char(x),'builtin') && ~exist(char(x),'file');