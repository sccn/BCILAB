function res = is_undefined_function(exp)
% test if a given function handle has no associated code (i.e. it is neither a lambda nor a reference to a known m-file)
res = false;
cexp = char(exp);
if cexp(1) ~= '@' && ~exist(cexp,'builtin')
    try
        % this section is a performance optimization due to the lack of a function exist(...,'mfileinpath'); the code is correct in all cases.
        % * for undefined functions, this line yields immediately a MATLAB:UndefinedFunction exception
        % * for all other functions, this line either yields, in almost all cases (excluding varargin+vararout functions),
        %   immediately a 'MATLAB:TooManyOutputs' or 'MATLAB:TooManyInputs' (and the like) exception
        %   or in the remaining cases either some exception further into the function's code or it succeeds.
        z=@someuncommontype;
        [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,t,u,v,w,x,y] = exp(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z,z); %#ok<NASGU,ASGLU>
    catch e
        % side note: below we make sure that the UndefinedFunction error is generated here and not deeper into a potentially defined exp()'s code
        res = strcmp(e.identifier,'MATLAB:UndefinedFunction') && strcmp(e.stack(1).name,mfilename);
    end
end