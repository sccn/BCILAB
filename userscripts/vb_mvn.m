function [mu_raw,prec_raw] = vb_mvn(varargin)
% Approximate the given posterior distribution with a multivariate Gaussian.

% test code:
% vb_mvn({'mu',3*ones(1,10),'sig',eye(12),'beta',3,'gamma',zeros(3)}, {@(mu,sig,beta) mu*sig+beta,@(mu,gamma)mu-gamma});
% vb_mvn({'beta',zeros(F+1,1)},{@(beta)logmvnpdf(y,Xm*beta(2:end)+beta(1),s2y*eye(T)),@(beta)logmvnpdf(beta,zeros(F+1,1),s2beta*eye(F+1))});

% parse arguments
opts = arg_define(varargin,...
    arg_norep({'variables','Variables'},[],[],'Variables of the problem. This is a struct or cell array of NVPs of the variables over which a posterior distribution shall be inferred. Each variable must be given with an initial value.','type','expression'),...
    arg_norep({'factors','Factors'},[],[],'Factors of the log-posterior. This is a cell array of anonymous functions, each of which represents one factor in the distribution. This is mostly for convenience to make it easier to specify the posterior. The arguments of the anonymous functions that should be inferred over should match the variable names.','type','expression'), ...
    arg({'vb_iters','NumVBIterations'},100,[],'Number of iterations for variational loop.'));
[variables,factors,vb_iters] = arg_toworkspace(opts);

% --- parse the variables of the problem ---

% get cell array of names and initial values for each variable
if isstruct(variables)
    var_names = fieldnames(variables)';
    var_init = struct2cell(variables);
elseif iscell(variables) && iscellstr(variables(1:2:end)) && mod(length(variables),2)==0
    var_names = variables(1:2:end);
    var_init = variables(2:2:end);
end

% a map from variable name to index into various arrays
varname2idx = cell2struct(num2cell(1:length(var_names)),var_names,2);

% get the number of elements for each variable and their shapes
var_numels = cellfun('prodofsize',var_init);
var_shapes = cellfun(@size,var_init,'UniformOutput',false);

nVars = length(var_names);           % number of formal variables
nDims = sum(var_numels);             % dimensionality of the parameter space
nHessNonzeros = sum(var_numels.^2);  % max. number of nonzeros in the log-posterior's Hessian

% get the index ranges that each variable occupies in the posterior parameter space
var_ranges = {}; [var_ranges{1:nVars}] = chopdeal(1:nDims,var_numels);

% --- parse the factors of the posterior ---

% determine which variables each factor depends on...
for f=length(factors):-1:1 %#ok<FXUP>
    % get the names of the factor's arguments
    str = char(factors{f});
    fact_names{f} = hlp_split(str(find(str=='(',1)+1:find(str==')',1)-1),',');
    fact_nArgs(f) = length(fact_names{f});
    % get remaining properties
    nArgs = fact_nArgs(f);
    for n=nArgs:-1:1
        varname = fact_names{f}{n};
        varidx = varname2idx.(varname);
        fact_numels{f}(n) = var_numels(varidx); % number of elements of the n'th argument of factor f
        fact_masks{f}{n} = var_ranges{varidx};  % index range of the argument in the full posterior parameter space
        fact_shapes{f}{n} = var_shapes{varidx}; % shape of the argument
    end
    fact_mask{f} = [fact_masks{f}{:}];             % index range to extract all arguments from the full posterior
    % assemble the Hessian mask for the factor from blocks of index ranges
    hess_masks{f} = cell(1,nArgs^2);
    % append diagonal block indices...
    for n=1:nArgs        
        tmp = spalloc(nDims,nDims,fact_numels{f}(n)^2);
        tmp(fact_masks{f}{n},fact_masks{f}{n}) = ones(fact_numels{f}(n)); %#ok<SPRIX>
        hess_masks{f}{n} = find(tmp);
    end
    % append off-diagonal block indices...
    k = n+1;
    for i=1:nArgs
        for j=[1:i-1,i+1:nArgs]
            tmp = spalloc(nDims,nDims,fact_numels{f}(i)*fact_numels{f}(j));
            tmp(fact_masks{f}{i},fact_masks{f}{j}) = ones(fact_numels{f}(i),fact_numels{f}(j)); %#ok<SPRIX>
            hess_masks{f}{k} = find(tmp);
            k = k+1;
        end
    end
    hess_mask{f} = vertcat(hess_masks{f}{:});    
end

% generate code for derivative calculation
for f=length(factors):-1:1 %#ok<FXUP>
    func = functions(factors{f});
    code_signature = hlp_cryptohash(func.function);    
    size_signature = num2str(hlp_fingerprint({fact_shapes{f},cellfun(@size,struct2cell([func.workspace{:}]),'UniformOutput',false)}));
    src_file = env_translatepath(['bcilab:/temp/vb/anon_' code_signature '_s' size_signature '.m']);
    if ~exist(src_file,'file')
        io_mkdirs(src_file,{'+w','a'});    
        try
            fid = fopen(src_file,'w');
            if fid == -1
                error('Cannot open file %s for writing.',src_file); end            
        catch e
        end
    end
end

% stack all initial values to get the initial posterior mean and precision matrix
mu_init = vertcat_vec(var_init);
prec_init = speye(length(mu_init));

% TODO: optionally do gradient descent and use the laplace approximation to initialize
%       mu_init and prec_init; I suspect that this will make the convergence of the rest
%       a lot more efficient

% do inference
[mu_raw,prec_raw] = GaussVarApproxHessian(mu_init,prec_init,@gradFun,vb_iters);

% TODO: build output variables



    function [grad,hess] = gradFun(beta)
        % evaluate the gradient and Hessian of the log posterior
        % [Gradient,Hessian] = gradFun(Beta)
        %
        % In:
        %   Beta : parameter setting for which to evaluate the log-posterior
        %
        % Out:
        %   Gradient : gradient of the posterior at Beta
        %
        %   Hessian : Hessian matrix of the posterior at Beta
        
        % initialize gradient and Hessian to zero
        grad = zeros(nDims,1);
        hess = spalloc(nDims,nDims,nHessNonzeros);
        
        % for each factor...
        for f=1:length(factors) %#ok<FXUP>
            nArgs = fact_nArgs(f);
            % read out the parameter subset from beta needed by factor f and reshape it
            X = cell(1,nArgs); 
            [X{:}] = chopdeal(beta(fact_mask{f}),fact_numels{f});
            for a=nArgs:-1:1
                X{a} = reshape(X{a},fact_shapes{f}{a}); end
            % evaluate factor
            if nargout == 1
                % ... only gradients
                [dummy,gradients{1:nArgs}] = factors{f}(X{:}); %#ok<ASGLU>
            else
                % ... including all Hessian sub-blocks (first the diagonals, then off-diagonals in reading order)
                [dummy,gradients{1:nArgs},hessians{1:nArgs^2}] = factors{f}(X{:}); %#ok<ASGLU>
            end
            % add to gradient and optionally Hessian to
            grad(fact_mask{f}) = grad(fact_mask{f}) + vertcat_vec(gradients);
            if nargout > 1
                hess(hess_mask{f}) = hess(hess_mask{f}) + vertcat_vec(hessians); end
        end        
    end


end

% vectorize and vertically concatenate the values in a given cell array
function v = vertcat_vec(v)
    v = cellfun(@vec,v,'UniformOutput',false);
    v = vertcat(v{:});
end
