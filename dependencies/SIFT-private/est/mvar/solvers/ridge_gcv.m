function [x,lambda_opt,svd_state,Y_hat] = ridge_gcv(varargin)
%[x_hat,lambda_opt] = ridgeGCV(Y,A,P)
%
% Estimates a ridge regression model, also know as Tikhonov regularization,
% or minimum norm with L2 prior.
%
% x_hat = argmin(x) ||Y-A*x||^2 + lambda*||P*x||^2
% with lambda > 0
%
% 
% [..., svd_state] = ridge_gcv(...)  returns the svd decomposition. This can
% be reused if the target vector does not change.
%
% [..., Y_hat] = ridge_gcv(...)  returns the predicted target vector.
% Residuals can be computed via E = Y-Y_hat;
%
% This code is based on a previous implementation used in Valdes-Hernandez
% et al. (2009), written by Alejandro Ojeda and Pedro Valdez-Hernandez at
% the Cuban Neuroscience Center in 2009.
%
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jul-2012
%         Tim Mullen, SCCN/INC/UCSD, Jan-2013, Apr-2013
%
% References:
%   Pedro A. Valdes-Hernandez, Alejandro Ojeda, Eduardo Martinez-Montes, Agustin
%       Lage-Castellanos, Trinidad Virues-Alba, Lourdes Valdes-Urrutia, Pedro A.
%       Valdes-Sosa, 2009. White matter architecture rather than
%       cortical surface area correlates with the EEG alpha rhythm. NeuroImage 49
%       (2010) 2328–2339

arg_define([0 Inf],varargin, ...
    arg_norep({'Y','TargetVector','y'},mandatory,[],'The target vector'), ...
    arg_norep({'A','DesignMatrix'},mandatory,[],'The design matrix. This is the data matrix (ie X).'), ...
    arg({'P','PriorInvCov'},[],[],'Prior inverse covariance (precision) matrix for params. Can be an [nc x nc] matrix, where nc is the number of columns of A. Can also be a scalar, P, specifying the prior inverse variance (precision) of each parameter (diagonal covariance matrix). If empty, identity covariance matrix assumed. A sparse matrix is advised if precision matrix is not dense.'), ...
    arg_nogui({'blksz','DesignMatrixBlockSize','designMatrixBlockSize'},[],[],'Design matrix structure. Can be a tuple [numrows numcols], in which case A consists of identical blocks of this size, along the main diagonal'), ...
    arg_subswitch({'lambdaMode','LambdaSelectionMode'},'grid_gcv', { ...
        'manual' { ...
            arg({'lambda','RegularizationParam'},1,[0 Inf],'Regularization parameter (lambda)','type','denserealdouble') ...
            }, ...
         'grid_gcv' { ...
            arg({'gridSize','GridSize'},100,[0 Inf],'Grid size for regularization param search. This is used to automatically select the regularization parameter which minimizes the Generalized Cross-Validation (GCV) criteria.') ...
            arg({'plotGCV','PlotGCV'},false,[],'Plot GCV curve'), ...
            } ...
     },'Selection mode for lambda. Automatic (GCV grid search) or Manual (must provide lambda)'), ...
    arg_nogui('svd_state',[],[],'SVD decomposition structure'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output') ...
    );


[nr,nc] = size(A);

% if A is not sufficiently sparse, convert to full
if issparse(A) && nnz(A)/numel(A) > 0.9
    A = full(A);
end

if isempty(svd_state)
    if verb, fprintf('Computing SVD of design matrix.\nr'); end
    
    % init prior covmat
    if isscalar(P)
        P=P*speye(nc);
    end
    if ~isempty(P)
        Pinv = inverse(P);
    end
    % compute SVD
    if ~isempty(P)
        APinv = A*Pinv;
    else
        APinv = A;
    end
    if isempty(blksz)
        [U,S,V] = svd_wrapper(APinv);
    else
        %warning('block optimization not yet implemented');
        [U,S,V] = svd_wrapper(APinv);
%         [U,S,V] = svd_wrapper(APinv(1:blksz(1),1:blksz(2)));
    end
    if isempty(P)
        iPV = V;
    else
        iPV = Pinv*V;
    end
    s   = diag(S);
    s2  = s.^2;
    Ut  = U';
else
    iPV = svd_state.iPV;
    s   = svd_state.s;
    s2  = svd_state.s2;
    Ut  = svd_state.Ut;
end

UtY = Ut*Y;

switch lambdaMode.arg_selection
    case 'grid_gcv'
        % search over a grid of lambda values for the value that minimizes
        % the Generalized Cross Validation (GCV) criteria
        % lambdaOpt = argmin(lambda) { GCV(lambda) }
        
        % automatically determine lambda range based on singular values
        tol     = max([nr nc])*eps(max(s));
        lgrid   = logspace(log10(tol),log10(max(s)),lambdaMode.gridSize);
        gcv     = zeros(lambdaMode.gridSize,1);
        for it=1:lambdaMode.gridSize
            % compute GCV criteria
            d       = lgrid(it)./(s2+lgrid(it));
            f       = diag(d)*UtY;
            gcv(it) = dot(f,f,1)/sum(d)^2;
        end
        loc = getMinima(gcv);
        if isempty(loc),
            % no minimum found, search for elbow instead
            if verb
                fprintf('no GCV minimum, finding elbow...\nr'); 
            end
            [val loc] = hlp_findElbow(gcv);  % min(gcv)
        end
        loc         = loc(end);
        lambda_opt  = lgrid(loc);
        if verb
            fprintf('lambda: %0.5g, GCV: %0.5g\nr',lambda_opt,gcv(loc)); 
        end
        if lambdaMode.plotGCV
            plotLambdaGCV(lgrid,gcv,lambda_opt,loc); 
        end
    case 'manual'
        lambda_opt = lambdaMode.lambda;
    otherwise
        error('SIFT:ridge_gcv:BadLambdaRule', ...
              'Unknown lambda learning rule %s',lambdaMode.arg_selection);
end

% solve system for parameters x
x = iPV*bsxfun(@times,(s./(s2+lambda_opt^2)),UtY);

if nargout > 2
    Y_hat = A*x;
end
if nargout > 3
    svd_state.iPV = iPV;
    svd_state.s   = s;
    svd_state.s2  = s2;
    svd_state.UtY = UtY;
end

% Helper functions
% -------------------------------------------------------------------------
function [U,S,V] = svd_wrapper(A)
% compute SVD
if issparse(A)
    [U,S,V] = svds(A,min(size(A)));
else
    [U,S,V] = svd(A,'econ');
end
        
function plotLambdaGCV(lgrid,gcv,lambda_opt,loc)
% plot lambda versus GCV
figure;
semilogx(lgrid,gcv)
xlabel('log-lambda');
ylabel('GCV');
hold on;
plot(lambda_opt,gcv(loc),'rx','linewidth',2);
hold off; grid on;


function indmin = getMinima(x)
% get minimum of a function
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
