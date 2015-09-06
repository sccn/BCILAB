function MCMC_State = stat_ubsplsm_init(varargin)
% Initialize MCMC estimator for univariate B-Spline smoother
%
% Author: Tim Mullen and Wes Thompson, 2010-13, SCCN/INC, UCSD.
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

arg_define([0 Inf],varargin, ...
    arg_norep({'Y','TSData'},mandatory,[],sprintf(['Cell array of data to smooth.\n' ...
              'Generally, Y is a cell array where Y{i} is the T x 1 vector of time-varying (or freq-varying) connectivity for the kth channel pair of the sth subject.\n' ...
              'e.g. Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ... \n' ...
              '          s2(1,1) s2(1,2) s2(1,3) ... } \n'])), ...
    arg({'K','Knots'},5,[],'Positions of spline knots. If K is a scalar, then we automatically assign K evenly spaced knots. A good heuristic is one knot every 5%','shape','row'), ...
    arg({'Q','FPCABasisDim','fpcaBasisDim'},4,[0 Inf],'Number of FPCA basis functions.'), ...
    arg({'order','BSplineOrder'},4,[1 Inf],'B-Spline order. 4 is a good choice'), ...
    arg({'initNoiseVariance'},0.1,[eps Inf],'Initial noise variance'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
);
% arg({'T','TimeSeriesLength'},mandatory,[],'Time series length. This is the length of the time-series (vector) to smooth'), ...
% arg({'R','NumTimeSeries'},mandatory,[],'Number of observations. This is the total number of time-series to smooth. For instance, this might be the total number of connectivity edges across all subjects.'), ...

T = size(Y{1},1);
R = length(Y);

% Determine knot locations
% -------------------------------------------------------------------------
if isscalar(K)
    idx   = 1:T;
    K = K+order-1;
    knots = quantile(idx,(0:1:(K-order+1)) / (K-order+1));
    knots = knots(:)';
else
    % knots are already provided
    knots = K(:);
    K = length(knots)+2;  % add endpoints
end
   
% Construct orthonormal basis functions
% -------------------------------------------------------------------------
phi_t = stat_ubsplsm_mkspl(knots,T,order,verb);

% Initialize FPCA
% -------------------------------------------------------------------------
if verb==2
    multiWaitbar('Initializing FPCA','Reset','Color',hlp_getNextUniqueColor);
end

% Initialize the FPCA components to random, orthonormal vectors
% [!] should be able to replace this whole section with a single line of orth(mvrnd(...)). We only need orthonormal random vectors
THETA(:,1,1) = mvnrnd(zeros(K,1),eye(K));
THETA(:,1,1) = THETA(:,1,1)/norm(THETA(:,1,1));
for q = 2:Q                                                                     
    if verb==2
        multiWaitbar('Initializing FPCA',q/Q);
    end
    THETA(:,q,1) = mvnrnd(zeros(K,1),eye(K));

    % orthonormalization of random vector
    for r = 1:(q-1)

        THETA(:,q,1) = THETA(:,q,1) ...
                       -(THETA(:,q,1)'*THETA(:,q-r,1)/norm(THETA(:,q-r,1).^2)) ...
                       * THETA(:,q-r,1);
    end
    THETA(:,q,1) = THETA(:,q,1)/norm(THETA(:,q,1));
end
    
% Construct MCMC state object
% -------------------------------------------------------------------------
MCMC_State.THETA     = THETA;
MCMC_State.ALPHA     = zeros(Q,R);
MCMC_State.ALPHA_BAR = zeros(Q,1);
MCMC_State.SIGMA_EPS = initNoiseVariance;
MCMC_State.phi_t     = phi_t;
MCMC_State.initstate = true;

if verb==2
    multiWaitbar('Initializing FPCA', 'Close');
end