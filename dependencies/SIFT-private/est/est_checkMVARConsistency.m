function stats = est_checkMVARConsistency(varargin)
%
% For a VAR[p] process fit to N time windows, this function returns a
% vector of MVAR model consistency estimates for each time window. The
% percent consistency [1] is an index of the ability for a VAR model, fit
% to data X, to generate data with the same covariance structure as X.
% If Rr and Rs represent the vectorized corss-correlation matrices of the
% real and simulated data, respectively, then the percent consistency is
% given by:
%
% PC = (1 - norm(Rs - Rr)/norm(Rr)) * 100
%
% A value near 100% indicates the model has well-captured the covariance
% structure of the original data
%
% Inputs:
%
%   EEG:            EEG dataset
%   MODEL:          MODEL structure
%   typeproc:       reserved for future use, set to 0 for now
%
% Outputs:
%
%   stats
%       .PC:             Percent Consistency Index [1]
%
% References:
%
% [1] Ding M, Bressler SL, Yang W, Liang H (2000) Short-window spectral
% analysis of cortical event-related potentials by adaptive multivariate
% autoregressive modeling: data preprocessing, model validation, and
% variability assessment. Biol. Cybern. 83:35-45
%
% [2] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
% Theoretical Handbook and User Manual. Chapters 3,6.
% Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% See Also: pop_est_validateMVAR(), est_consistency()
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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

g = arg_define([0 2],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_norep({'MODEL','Model'},mandatory,[],'MVAR MODEL object'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model. Default is empty (use all windows)','cat','Options'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Options'), ...
    arg({'Nr','NumRealizations'},[],[],'Number of realizations (trials) to simulate. Default is the larger of 30 and the number of trials in the dataset.','type','int32'), ...
    arg({'donorm','Normalize'},false,[],'Normalize data. This will z-score (standardize) both simulated and real data before computing correlations'), ...
    arg({'nlags','NumLags'},[],[],'Number of correlation lags. If empty, max(10,model_order) is used.','type','int32'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

% commit EEG and MODEL variables to workspace
[data g] = hlp_splitstruct(g,{'EEG','MODEL'});
arg_toworkspace(data);
clear data;

% determine correlation lags
if isempty(g.nlags)
    g.nlags = max(10,MODEL.morder); 
end
minLagAcf = MODEL.morder;
maxLagAcf = floor(MODEL.winlen*EEG.srate)-MODEL.morder-1;
if (g.nlags < minLagAcf)
    fprintf('checkMVARConsistency: Number of correlation lags is less than mininum allowed (%d lags). Will use %d lags for checks.\n',minLagAcf,minLagAcf);
    g.nlags = minLagAcf;
elseif (g.nlags > maxLagAcf)
    fprintf('checkMVARConsistency: Number of correlation lags exceeds maximum allowed (%d lags). Will use %d lags for checks.\n',maxLagAcf,maxLagAcf);
    g.nlags = maxLagAcf;
end

% window size in points
winLenPnts = floor(MODEL.winlen*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round(MODEL.winStartTimes*EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    g.winArrayIndex = randwin;
end

% get the array indices of the windows we are working with
g.winArrayIndex = getindex(MODEL.winStartTimes,(g.winStartIdx-1)/EEG.srate);

% initialize waitbar
if g.verb==2
    waitbarTitle = sprintf('Checking consistency %s...', ...
        fastif(isempty(EEG.condition),'',['for ' EEG.condition]));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle, ...
                 'Color', hlp_getNextUniqueColor, ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
%     wb_cleanup = onCleanup(@() multiWaitbar(waitbarTitle,'Close'));
end

numWins = length(g.winStartIdx);
stats.PC = zeros(1,numWins);

for t=1:numWins
    
    % get the array index of the window we are working with
    winArrIdx = g.winArrayIndex(t);
    if isfield(MODEL,'mu') && ~isempty(MODEL.mu)
        mu = MODEL.mu{winArrIdx};
    else
        mu = [];
    end
    data = squeeze(EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:));
    stats.PC(t)= est_consistency(data,MODEL.AR{winArrIdx},MODEL.PE{winArrIdx},g.donorm,g.Nr,g.nlags,mu);
    
    if g.verb==2
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,t/numWins);
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            stats = [];
            return;
        end
    end
    
end

stats.winStartIdx   = g.winStartIdx;
stats.winStartTimes = MODEL.winStartTimes(g.winArrayIndex);
stats.winArrayIndex = g.winArrayIndex;

% clean up
if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end


