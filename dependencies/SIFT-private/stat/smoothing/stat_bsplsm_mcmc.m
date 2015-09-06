

function GAMMA_HAT = stat_bsplsm_mcmc(C, bsplmodel,niter,zeta,pi1,pi2)

% FUNCTION UNDER CONSTRUCTION

% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
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

% C 
%
% bsplmodel is the b-spline basis model from stat_bsplsm_mkspl
%
% niter is the number of MCMC iterations
%
% niterToKeep is the number of (final) iterations to return in GAMMA_HAT
% 
% g = arg_define([0 1],varargin, ...
%     arg_norep({'C','TSData'},mandatory,[],'Connectivity Matrix. This is a [freqs x times] matrix of connectivity values for a pair of variables/channels'), ...
%     arg_nogui({'bsplmodel','BSplineModel'},struct([]),[],'Optional bivariate B-Spline model. This can be computed via stat_bsplsm_mkspl().'), ...
%     arg({'K','Knots'},5,[],'Positions of spline knots along frequency dimension (Hz). If K is a scalar, then K knots are evenly spaced from first to last frequency. A good heuristic is one knot every 5%','shape','row'), ...
%     arg({'Q','FPCABasisDim','fpcaBasisDim'},4,[0 Inf],'Number of FPCA basis functions.'), ...
%     arg({'smoothingLayout','MatrixElementsToSmooth'},{'diagonals','off-diagonals'},{'off-diagonals'},'Which parts of the matrix to smooth. Diagonals (e.g. auto-connectivity) and off-diagonals (e.g. cross-connectivity) will be smoothed separately','type','logical'), ...
%     arg({'niters','nMCMCiters','NumMcmcIters','niter'},1000, [1 Inf], 'Number of MCMC iterations for spline fitting'), ...
%     arg({'basisCoeffVarPrior'},1000,[eps Inf],'Variance of basis coefficient gaussian prior. Larger --> more wiggling allowed'), ...
%     arg({'noiseVarPriorShape'},0.01,[eps Inf],'Shape (D.O.F) of noise variance prior. This is the "alpha" parameter of the inverse gamma prior distribution. Increasing noiseVarPriorShape --> decreased variance of noise variance distribution.'), ...
%     arg({'noiseVarPriorScale'},0.01,[eps Inf],'Scale parameter of noise variance prior. This is the "theta" (1/beta) parameter of inverse gamma prior distribution. Increasing noiseVarPriorScale --> right-shift of distribution --> (increase in expected noise variance). In general MEAN(noiseVariance) = noiseVarPriorScale/noiseVarPriorShape and MODE(noiseVariance) = noiseVarPriorScale/(noiseVarPriorShape-1) for noiseVarPriorShape>=1.'), ...
%     arg({'initNoiseVariance'},0.1,[eps Inf],'Initial noise variance'), ...
%     arg_norep({'MCMC_InitState'},struct([]),[],'Object containing initial state of Gibbs sampler'));
%     arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
% 



if nargin<3
    niter = 1000;
end

% commit bspline model to workspace 
arg_toworkspace(bsplmodel);

%% Set up hyperparameters

% prior 
if ~exist('zeta','var')
    zeta=1000; end

% prior
if ~exist('pi1','var')
    pi1=1; end

% prior
if ~exist('pi2','var')
    pi2=1; end


%% Initialize parameters
N=size(Phi_R,1);
NPhiF=size(Phi_F,1);
Q = 4;

H=H1*H2;
ETA=zeros(Q,niter+1);
DELTA=zeros(H-Q,niter+1);
SIGMA_SQ_EPS=zeros(1,niter+1);
SIGMA_SQ_DELTA=zeros(2,niter+1);
GAMMA_HAT=zeros(ntimes,nfreqs,niter+1);

C = C';   % make C be [times x freqs]

% Starting values

ETA(:,1)     = normrnd(0,0.1,[Q 1]);
DELTA(:,1)   = normrnd(0,0.1,[H-Q 1]);

SIGMA_SQ_EPS(1) = 1;
SIGMA_SQ_DELTA(:,1)  = ones(2,1);


for iter=1:niter
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw ETA                                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inv_Sigma_eta=eye(Q)/zeta;
    INV_SIGMA_SQ_EPS = 1/SIGMA_SQ_EPS(iter);
    mu_eta=sum(reshape(Phi_F,NPhiF,ntimes*nfreqs).*repmat((C(:)'-DELTA(:,iter)'*reshape(Phi_R,N,ntimes*nfreqs))*INV_SIGMA_SQ_EPS,NPhiF,1),2);
    
    for t=1:ntimes
        for f=1:nfreqs
            inv_Sigma_eta=inv_Sigma_eta+squeeze(Phi_F(:,t,f))*squeeze(Phi_F(:,t,f))'*INV_SIGMA_SQ_EPS;
        end
    end
    
    Sigma_eta=double(inverse(inv_Sigma_eta));
    Sigma_eta=.5*(Sigma_eta+Sigma_eta');
    mu_eta=Sigma_eta*mu_eta;
    ETA(:,iter+1)=mvnrnd(mu_eta,Sigma_eta)';
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw DELTA                                 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inv_Sigma_delta=SIGMA_SQ_DELTA(1,iter)\S1_star+SIGMA_SQ_DELTA(2,iter)\S2_star;
    
    mu_delta=sum(reshape(Phi_R,N,ntimes*nfreqs).*repmat((C(:)'-ETA(:,iter+1)'*reshape(Phi_F,NPhiF,ntimes*nfreqs))*INV_SIGMA_SQ_EPS,N,1),2);
    
    for t=1:ntimes 
        for f=1:nfreqs
            inv_Sigma_delta=inv_Sigma_delta+squeeze(Phi_R(:,t,f))*squeeze(Phi_R(:,t,f))'*INV_SIGMA_SQ_EPS;
        end
    end
    Sigma_delta= double(inverse(inv_Sigma_delta));
    Sigma_delta=.5*(Sigma_delta+Sigma_delta');
    mu_delta=Sigma_delta*mu_delta;
    DELTA(:,iter+1)=mvnrnd(mu_delta,Sigma_delta)';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               compute GAMMA_HAT                          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    GAMMA_HAT(:,:,iter+1) = reshape(ETA(:,iter+1)'*reshape(Phi_F,NPhiF,ntimes*nfreqs)+DELTA(:,iter+1)'*reshape(Phi_R,N,ntimes*nfreqs),ntimes,nfreqs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw SIGMA_SQ_EPS                          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pi1_eps = pi1+.5*nfreqs*ntimes;
    pi2_eps = 0.5*(C-GAMMA_HAT(:,:,iter+1)).^2;
    pi2_eps = sum(pi2_eps(:))+pi2;
    
    
    SIGMA_SQ_EPS(iter+1)=gamrnd(pi1_eps,1/pi2_eps);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               draw SIGMA_SQ_DELTA                        %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pi1_delta=pi1+.5*(H-Q);
    pi2_delta=pi2+.5*DELTA(:,iter+1)'*(S1_star+S2_star)*DELTA(:,iter+1);
    SIGMA_SQ_DELTA(:,iter+1)=repmat(gamrnd(pi1_delta,1/pi2_delta),1,2);
    
    
end


% retain only final nitersToKeep of MCMC iterations
nd = ndims(GAMMA_HAT);    
switch nd
    case 3
        GAMMA_HAT = GAMMA_HAT(:,:,end-nitersToKeep+1:end);
    case 2
        GAMMA_HAT = GAMMA_HAT(:,end-nitersToKeep+1:end);
end
            


%%

% smoothY = mean(GAMMA_HAT(:,:,round(iter/2):end),3);
% 
% figure;
% subplot(1,2,1),surf(times,freqs,C'); xlabel('Time'); ylabel('Freq');
% subplot(1,2,2),surf(times,freqs,smoothY'); xlabel('Time'); ylabel('Freq');
% 
% 
% baseline = [1 10];
% basedist = squeeze(mean(GAMMA_HAT(:,baseline,round(iter/2):end),2));
% thresh = prctile(basedist,95,2);
% smoothY(smoothY<repmat(thresh,1,size(smoothY,2))) = 0;
% figure;
% imagesc(times,freqs,smoothY'); set(gca,'Ydir','normal'); xlabel('Time'); ylabel('Freq');
% 


