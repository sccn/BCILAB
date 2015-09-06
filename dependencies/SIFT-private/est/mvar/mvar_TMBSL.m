function [AR PE lambdaOpt] = mvar_TMBSL(varargin)
% Algorithm: T-MSBL
%
% Description:
%
% n/a
%
% Dependencies: TMSBL()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization coefficient
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix
%
% Requires: 
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), est_fitMVARKalman()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual.
%     Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen 2012, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

verb = arg_extract(varargin,{'verb','Verbosity'},[],0);

g = arg_define(varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'p','ModelOrder','morder'},10,[],'VAR Model order'), ...
                arg({'verb','Verbosity'},verb,[],'Verbose output') ...
                );

%             arg({'noise','NoiseLevel'},'custom',{'custom','no','small','mild','large'},'Noise Level'), ...
%                 arg({'PRUNE_GAMMA','PruneGamma'},1e-3,[0 Inf],sprintf(['Pruning threshold.\n' ...
%                                                                        'Threshold for prunning small hyperparameters gamma_i.\n'    ...
%                                                                        'In noisy cases, you can set MIN_GAMMA = 1e-3 or 1e-4. \n'   ...
%                                                                        'In strong noisy cases (e.g. SNR <= 6 dB), set MIN_GAMMA = 1e-2 for better performance.' ...
%                                                                        ])), ...
%                 arg({'lambda','InitialReguParam','ReguParam','Lambda'},1e-3,[0 Inf],'Initial regularization value (lambda).','type','denserealdouble'), ...
%                 arg({'Learn_Lambda','LearnLambda'},true,[],'Learn regularization param Lambda. If set, use the lambda learning rule (thus the input lambda is just the initial value). But note the learning rule is not very good if SNR <= 6 dB (but still has better performance than some other algorithms)'), ...
%                 arg({'Enhance_Lambda','EnhanceLambda'},true,[],'Enhance regularization param Lambda. In large/mild noisy cases (<=22dB), set EnhanceLambda for better performance.'), ...
%                 arg({'Matrix_Reg','Matrix_Reg'},1,[0 Inf],'Parameter covariance prior (sigma). If sigma a scalar, the parameter covariance prior is taken to be a diagonal matrix with sigma (variance) on diagonals. Otherwise, sigma can be a full prior covariance matrix. A sparse matrix is advised if covariance matrix is not dense.'), ...
%                 
            
% arg({'scaled','UseScaling'},true,[],'Scaling option. If set, coefficient estimates are restored to the scale of the original data'), ...

arg_toworkspace(g);

[nchs npnts ntr] = size(data);

% assemble predictors X and target variables Y for the structural equation
% Y = XA + noise

X = [];
Y = [];
for itr = 1:ntr
    Xi = [];
    for jj = 1:p
        Xi = cat(3, Xi, squeeze(data(:, p+1-jj:end-jj, itr))');
    end
    X = cat(1, X, reshape(permute(Xi, [1 3 2]), npnts-p, nchs*p));
    Y = cat(1,Y,data(:, p+1:end,itr)');
end

Y = vec(Y);
X = blkdiageye(X,nchs);

switch g.lambdaSelMode.arg_selection
    case 'manual'
        lambdaOpt = g.lambdaSelMode.lambda;
        gridSize  = 0;
    case 'automatic'
        gridSize  = g.lambdaSelMode.lambdaGridSize;
        lambdaOpt = [];
end
       
if isscalar(g.Matrix_Reg)
    g.Matrix_Reg=g.Matrix_Reg*speye(size(X,2));
end

[H2] = TMSBL(X,Y);

H2 = reshape(H2,[p,nchs,nchs]);

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res',1);
end


function v = vec(x)
v = x(:);
