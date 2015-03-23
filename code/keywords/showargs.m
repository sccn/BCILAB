function showargs(funcname)
% Show arguments of an argument-aware function with their descriptions
% Usage: showargs funcname

% query arguments
args = arg_report('lean',str2func(funcname));

% print report
fprintf('Arguments of funcname: \n');
for a=args    
    if length(a.names) > 1
        name = a.names{2};
    else
        name = a.names{1};
    end
    fprintf(' %s -- %s (%s,%s)\n',name, a.help{1}, a.type, a.shape);
end
