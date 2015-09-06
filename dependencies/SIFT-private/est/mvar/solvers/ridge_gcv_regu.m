function [x,lambda_opt,svd_state,Y_hat] = ridge_gcv_regu(varargin)
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
% Author: Tim Mullen, SCCN/INC/UCSD, Jan-2013, Apr-2013
%         Uses routines from Per Christian Hansen's Regularization Tools
%         package [1-2].
%
% References:
% [1] P. C. Hansen, Regularization Tools: A Matlab package for analysis and solution of discrete ill-posed problems, Numerical Algorithms, 6 (1994), pp. 1-35.
% [2] http://www.imm.dtu.dk/~pcha/Regutools/

arg_define([0 Inf],varargin, ...
    arg_norep({'Y','TargetVector','y'},mandatory,[],'The target vector'), ...
    arg_norep({'A','DesignMatrix'},mandatory,[],'The design matrix. This is the data matrix (ie X).'), ...
    arg_nogui({'blksz','DesignMatrixBlockSize','designMatrixBlockSize'},[],[],'Design matrix structure. Can be a tuple [numrows numcols], in which case A consists of identical blocks of this size, along the main diagonal'), ...
    arg_subswitch({'lambdaMode','LambdaSelectionMode'},'grid_gcv', { ...
        'manual' { ...
            arg({'lambda','RegularizationParam'},1,[0 Inf],'Regularization parameter (lambda)','type','denserealdouble') ...
            }, ...
         'grid_gcv' { ...
            arg({'plotGCV','PlotGCV'},false,[],'Plot GCV curve'), ...
            } ...
     },'Selection mode for lambda. Automatic (GCV grid search) or Manual (must provide lambda)'), ...
    arg_nogui('svd_state',[],[],'SVD decomposition structure'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output') ...
    );

% if A is not sufficiently sparse, convert to full
if issparse(A) && nnz(A)/numel(A) > 0.9
    A = full(A);
end

if isempty(svd_state)
    if verb, fprintf('Computing SVD of design matrix.\nr'); end
    
    % compute SVD
    if isempty(blksz)
        [U,s,V] = svd_wrapper(A);
    else
        %warning('block optimization not yet implemented');
        [U,s,V] = svd_wrapper(A);
%         [U,S,V] = svd_wrapper(APinv(1:blksz(1),1:blksz(2)));
    end
    s = diag(s);
else
    V = svd_state.V;
    s = svd_state.s;
    U = svd_state.U;
end

switch lambdaMode.arg_selection
    case 'grid_gcv'
        % search over a grid of lambda values for the value that minimizes
        % the Generalized Cross Validation (GCV) criteria
        % lambdaOpt = argmin(lambda) { GCV(lambda) }
        
        [lambda_opt,G,reg_param,minG] = gcv(U,s,Y,'Tikh');
        if verb
            fprintf('lambda: %0.5g, GCV: %0.5g\nr',lambda_opt,minG); 
        end
        if lambdaMode.plotGCV
            plotLambdaGCV(reg_param,G,lambda_opt,minG); 
        end
    case 'manual'
        lambda_opt = lambdaMode.lambda;
    otherwise
        error('SIFT:ridge_gcv_regu:BadLambdaRule', ...
              'Unknown lambda learning rule %s',lambdaMode.arg_selection);
end

% solve system for parameters x
x = tikhonov(U,s,V,Y,lambda_opt);

if nargout > 2
    Y_hat = A*x;
end
if nargout > 3
    svd_state.V = V;
    svd_state.s = s;
    svd_state.U = U;
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
        
function plotLambdaGCV(reg_param,G,reg_min,minG)
% plot lambda versus GCV
% Plot GCV function.
figure;
loglog(reg_param,G,'-'), xlabel('\lambda'), ylabel('G(\lambda)')
title('GCV function')
ax = axis;
HoldState = ishold; hold on;
loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
axis(ax)
if (~HoldState), hold off; end


function indmin = getMinima(x)
% get minimum of a function
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
