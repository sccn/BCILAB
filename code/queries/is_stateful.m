function r = is_stateful(f,varargin)
% check whether a filter function is stateful (for some arguments)
% stateful filter functions are characterized by a second output

nout = nargout(f);
if nout >= 0
    % usually, nargout tells us this
    r = nout > 1;
else
    % but in some cases (e.g., when the function is wrapped by a lambda), we must look more closely
    try
        [a,b] = f(varargin{:}); %#ok<NASGU,ASGLU>
        r = true;
    catch e
        r = ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs'}));
    end
end
