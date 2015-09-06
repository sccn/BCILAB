
function [stats residVar params] = est_checkMVARWhiteness(varargin)

% Tests for whiteness of residuals of VAR model. Residuals are 'white' if
% they are statistically uncorrelated (e.g. a white noise process). White
% residuals indicate the second-order dynamics of the data is adequately
% captured by the VAR model.
%
% If stats.pval > 0.05, we cannot reject the null hypothesis (that
% the data is white) at the 0.05 significance level, in other
% words the probability of mistakenly assuming the data is white is less
% than 5%, or, equivalently, there is at least a 95% probability that we 
% are correct in assuming the data is white.
%
% INPUT:
%
%   EEG:         EEG data structure
%   MODEL:       MODEL structure returned by est_fitMVAR or similar
%
% OPTIONAL:
%   'prct2sample':      {def: 100} percent of windows to sample [0 100]
%   'alpha':            {def: 0.05} whiteness significance level
%   'whitenessCriteria': Whiteness tests. Can be one or more of 
%                       'ACF','Ljung-Box','Box-Pierce','Li-McLeod'
%   'verb':             {def: true} enable verbosity
%   'numLagsAcf':       {def: 50} number of autocorrelation lags for 
%                       whiteness tests
% OUTPUT:
%   stats:
%    .[whitenessCriteria{i}]
%       .w:         vector of logicals. w(t)=1 if residuals for window t are 
%                   white, 0 otherwise
%       .pval:      vector of whitness p-values. If pval(t) < alpha,
%                   residuals are white at the alpha-level for window t.
%                   Alpha should be in [0 1]
%       .value:     values of the statistic
%       .fullname:  original name of test statistic
%       .winStartTime: start times of each window (sec)
%   racvf:     [nchs x nchs x 2*numLagsAcf+1 x numwins] matrix of autocorrelation
%            coeffs for lags 0:numLagsAcf
%   params:  struct of options used
%
%
% See Also: est_fitMVAR(), pop_est_validateMVAR()
%
% References:
%
% [1] Li, W. K., and A. I. McLeod, (1981). Distribution of the
%     Residual Autocorrelations in Multivariate ARMA Time Series
%     Models, J. Roy. Stat. Soc. B, 43, 231--239.      
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%     Springer.
% [3] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6. 
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
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
        arg({'alpha','SignificanceLevel'},0.05,[0 1],'Significance threshold. If a whitness test renders a (p-value > alpha), we cannot reject the null hypothesis (that the data is white) at the alpha-significance level, in other words the probability of mistakenly assuming the data is white is less than alpha, or, equivalently, there is at least a 1-alpha probability that we are correct in assuming the data is white.','cat','Statistics'), ...
        arg({'statcorrection','MultipleComparisonsCorrection'},'none',{'none'}, 'Correction for multiple comparisons. Not needed for multivariate tests','cat','Statistics'), ...
        arg({'numAcfLags','NumberOfAutocorrelationLags','NumACFLags'},50,[0 Inf],'Number of autocorrelation lags. This defines how many lags of the residual ACF must be jointly approx zero','cat','Statistics'), ...
        arg({'whitenessCriteria', 'WhitenessCriteria'}, {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'}, {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'},'Whiteness criteria. These are the statistical tests used to test for uncorrelated residuals','type','logical','cat','Statistics'), ...
        arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model. Default is empty (use all windows)','cat','Data Reduction'), ...
        arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Data Reduction'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
        );
       
% commit EEG and MODEL variables to workspace
[data g] = hlp_splitstruct(g,{'EEG','MODEL'});
arg_toworkspace(data);
clear data;

% do some input checking
minLagAcf = MODEL.morder;
maxLagAcf = floor(MODEL.winlen*EEG.srate)-MODEL.morder-1;
if (g.numAcfLags < minLagAcf)
    fprintf('WARNING: Number of autocorrelation lags is less than mininum allowed (%d lags). Will use %d lags for checks.\n',minLagAcf,minLagAcf);
    g.numAcfLags = minLagAcf;
elseif (g.numAcfLags > maxLagAcf)
    fprintf('WARNING: Number of autocorrelation lags exceeds maximum allowed (%d lags). Will use %d lags for checks.\n',maxLagAcf,maxLagAcf);
    g.numAcfLags = maxLagAcf;
end

% initialize outputs
[stats residVar] = deal([]);
if nargout > 2, params = g; end

if any(ismember_bc(lower(EEG.CAT.MODEL.algorithm),{'kalman','dekf'}))
    error('Whiteness tests currently not compatible with method ''%s''',EEG.CAT.MODEL.algorithm);
end

% convert whiteness criteria names into valid variable names
if ~isempty(g.whitenessCriteria)
    g.whitenessCriteria = lower(hlp_variableize(g.whitenessCriteria));
    if ~iscell(g.whitenessCriteria)
        g.whitenessCriteria = {g.whitenessCriteria}; 
    end
end
morder = MODEL.morder;

% convert window size to points
winLenPnts = round(MODEL.winlen*EEG.srate); 
[nchs npnts ntr] = size(EEG.CAT.srcdata);
% get total number of samples used to fit model
npntsTotal = ntr*(winLenPnts-morder);            

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

numWins = length(g.winStartIdx);

% initialize waitbar
if g.verb==2
    waitbarTitle = sprintf('%s %s...', ...
        fastif(isempty(g.whitenessCriteria),'Calculating residuals','Checking whiteness'), ...
        fastif(isempty(EEG.condition),'',['for ' EEG.condition]));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle, ...
                 'Color', hlp_getNextUniqueColor, ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end

for t=1:numWins
    
    % calculate residuals
    % ---------------------------------------------------------------------
    
    % get the array index of the window we are working with
    winArrIdx = g.winArrayIndex(t);
    
    % get the estimated mean of the process
    if ~isfield(MODEL,'mu') || isempty(MODEL.mu)
        mu = zeros(1,EEG.CAT.nbchan);
    else
        mu = MODEL.mu{winArrIdx};
    end
    
    % calc residuals
    residuals = est_mvarResiduals(squish(EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:)), ...
                                  MODEL.AR{winArrIdx}, mu,1,true);
    
                              
    % calculate residual autocovariance and autocorrelation matrices
    % ---------------------------------------------------------------------
    [racvf{t} racf{t}] = deal(zeros(ntr,nchs,nchs*(g.numAcfLags+1))); % residual autocov, autocorr
    resvar{t}          = zeros(ntr,nchs);  % residual variance
    
    % first calculate autocovariance/autocorrelation for each trial ...
    for tr=1:ntr
        % R = [Rs1s1 Rs1s2 Rs1s3 Rs2s1 Rs2s2 Rs2s3 Rs3s1 Rs3s2 Rs3s3 ...]
        tmpacvf = xcov(residuals(:,:,tr)',g.numAcfLags,'biased');

        % extract positive lags
        tmpacvf = tmpacvf(g.numAcfLags+1:end,:);         
        for lag = 1:g.numAcfLags+1
            racvf{t}(tr,:,(1:nchs)+nchs*(lag-1)) = reshape(tmpacvf(lag,:),[nchs nchs])';
        end

        % convert autocov to autocorrelation function
        rvar{t}(tr,:) = 1./sqrt(diag(squish(racvf{t}(tr,1:nchs,1:nchs))));
        D = diag(rvar{t}(tr,:));
        for lag = 1:g.numAcfLags+1
            racf{t}(tr,:,(1:nchs)+nchs*(lag-1)) = D*squish(racvf{t}(tr,:,(1:nchs)+nchs*(lag-1)))*D;
        end
    end
    % ... finally average autocovariance/correlation matrices across trials
    racf{t}      = squish(mean(racf{t},1));
    racvf{t}     = squish(mean(racvf{t},1));
    rvar{t}      = squish(mean(rvar{t},1));
    
    % compute whiteness criteria
    % ---------------------------------------------------------------------
    % TODO: performance can be improved by avoiding recomputation of similar statistics 
    for i=1:length(g.whitenessCriteria)
        switch lower(g.whitenessCriteria{i})
            
            case 'acf'
                % test residual autocorrelations
                sigthresh = 1.96/sqrt(npntsTotal);    % 0.95 confidence interval for acf
                numNonSigCoeffs = sum(sum(abs(racf{t}(:,nchs+1:end))<sigthresh));    % count how many coefficients are within the bounds for 'white' residuals
                numCoeffs = numel(racf{t}(:,nchs+1:end));                            % count the total number of coefficients
                stats.acf.pval(t) = numNonSigCoeffs / numCoeffs;                     % estimate the probability for a coefficient to be inside the bounds (probability of whiteness)
                stats.acf.w(t) = stats.acf.pval(t) > 1-g.alpha;                       
                stats.acf.acffun{t} = racf{t};
                stats.acf.fullname = 'ACF';
                stats.acf.winStartTimes = g.winStartIdx*EEG.srate;
            case 'ljungbox' 
                % ljung-box (modified portmanteau) test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(racvf{t}(1:nchs,1:nchs));
                for k=1:g.numAcfLags
                    Qh = Qh + 1/(npntsTotal-k) * trace(racvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*racvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.ljungbox.value(t) = Qh*npntsTotal*(npntsTotal+2);
                stats.ljungbox.pval(t) = 1-chi2cdf(stats.ljungbox.value(t),(nchs^2)*(g.numAcfLags-morder));
                stats.ljungbox.w(t)= stats.ljungbox.pval(t)>g.alpha;
                stats.ljungbox.fullname = 'Ljung-Box';
                stats.ljungbox.winStartTimes = g.winStartIdx*EEG.srate;
            case 'boxpierce' 
                % box-pierce portmanteau test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(racvf{t}(1:nchs,1:nchs));
                for k=1:g.numAcfLags
                    Qh = Qh + trace(racvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*racvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.boxpierce.value(t) = Qh*npntsTotal;   % Qh*npntsTotal^2  is the modified portmanteau test (Lutkepohl, p.171)
                stats.boxpierce.pval(t) = 1-chi2cdf(stats.boxpierce.value(t),(nchs^2)*(g.numAcfLags-morder));
                stats.boxpierce.w(t)= stats.boxpierce.pval(t)>g.alpha;
                stats.boxpierce.fullname = 'Box-Pierce';
                stats.boxpierce.winStartTimes = g.winStartIdx*EEG.srate;
            case 'limcleod'
                % Li-McLeod portmanteau test for residual autocorrelation up to lag h
                Qh=0;
                C0inv = inverse(racvf{t}(1:nchs,1:nchs));
                for k=1:g.numAcfLags
                    Qh = Qh + trace(racvf{t}(1:nchs,(1:nchs)+k*nchs)'*C0inv*racvf{t}(1:nchs,(1:nchs)+k*nchs)*C0inv);
                end
                stats.limcleod.value(t) = Qh*npntsTotal + nchs^2*g.numAcfLags*(g.numAcfLags+1)/(2*npntsTotal);
                stats.limcleod.pval(t) = 1-chi2cdf(stats.limcleod.value(t),(nchs^2)*(g.numAcfLags-morder));
                stats.limcleod.w(t) = stats.limcleod.pval(t) > g.alpha;    % using bonferroni correction for multiple tests (across channels)
                stats.limcleod.fullname = 'Li-McLeod';
                stats.limcleod.winStartTimes = g.winStartIdx*EEG.srate;
            otherwise
                fprintf('Unknown whiteness test %s',g.whitenessCriteria{i}); 
        end
    end
    
    if g.verb==2
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,t/numWins);
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            [stats, residVar] = deal([]);
            return;
        end
    end
                
end

% clean up
if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end

stats.winStartIdx        = g.winStartIdx;
stats.winStartTimes      = MODEL.winStartTimes(g.winArrayIndex);
stats.winArrayIndex      = g.winArrayIndex;
stats.alpha = g.alpha;

if nargout > 1
    residVar.autocorrelation = racf;
    residVar.autocovariance  = racvf;
    residVar.variance        = rvar;
    residVar.numAcfLags      = g.numAcfLags;
    residVar.winStartIdx     = g.winStartIdx;
    residVar.winStartTimes   = MODEL.winStartTimes(g.winArrayIndex);
    residVar.winArrayIndex   = g.winArrayIndex;
end
