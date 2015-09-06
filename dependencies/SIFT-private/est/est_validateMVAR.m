function [whitestats PCstats stabilitystats residualstats] = est_validateMVAR(varargin)
%
% Validate a fitted VAR model. With two inputs, this function generates a
% GUI where the validation scheme can be specified. Validation consists of
% statistical tests for "whiteness" of fitted VAR model residuals [1,2],
% consistency of the fitted model [1,3] and stability of fitted model [1-3]
%
% Input:
%
%   ALLEEG:     Array of EEGLAB data structures containing fitted MODEL
%   typeproc:   reserved for future use. Use 0
%
% Optional:
%
%     'whitenessCriteria':    Cell array containing names of residual whiteness
%                             test criteria to evaluate. See [1, 2] for details.
%                             Possible Values: {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'}
%                             Default Value  : all
%                             Data Input Type: cell array
%
%
%     'checkWhiteness':       Whether or not to check whiteness.
%                             Default Value  : true
%                             Data Input Type: boolean
%
%     'checkConsistency':     Whether or not to check consistency. See [1,3]
%                             for details.
%                             Default Value  : true
%                             Data Input Type: boolean
%
%     'checkStability':       Whether or not to check stability. See
%                             [1-3] for details.
%                             Default Value  : true
%                             Data Input Type: boolean
%
%     'alpha':                significance level for determining whiteness
%                             Data Input Range: [0 1]
%                             Default Value   : 0.05
%                             Data Input Type : real number (double)
%
%     'prctWinToSample':      percent of time windows to randomly select
%                             Data Input Range: [0 100]
%                             Default Value   : 100
%                             Data Input Type : real number (double)
%
%     'verb':                 verbosity level (0=no output, 1=text, 2=gui)
%
%
% Output:
%
%     whitestats:             Structure containing whiteness statistics.
%                             See est_checkMVARWhiteness() for details on
%                             structure format
%
%     PC:                     Vector of percent consistency estimates for
%                             each window.
%
%     stability:              Vector of stability estimates for each window
%
% See Also: est_checkMVARWhiteness(), est_checkMVARStability(),
%           est_checkMVARConsistency, pop_est_fitMVAR()
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 3.6 and 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%   Springer.
%
% [3] Ding M, Bressler SL, Yang W, Liang H (2000) Short-window spectral
%   analysis of cortical event-related potentials by adaptive multivariate
%   autoregressive modeling: data preprocessing, model validation, and
%   variability assessment. Biol. Cybern. 83:35-45
%
% Author: Tim Mullen, 2010, SCCN/INC, UCSD.
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

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],0);
prctWinToSample = arg_extract(varargin,{'prctWinToSample','WindowSamplePercent'},[],100);

g = arg_define([0 1],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_subtoggle({'checkWhiteness','CheckResidualWhiteness'},{'verb',verb,'prctWinToSample',prctWinToSample},@est_checkMVARWhiteness,sprintf('Check residual whiteness. \nResiduals are "white" if they are uncorrelated.'),'cat','Validation Methods','suppress',{'prctWinToSample','VerbosityLevel'}), ...
    arg_subtoggle({'checkResidualVariance','CheckResidualVariance'},{'verb',verb,'prctWinToSample',prctWinToSample,'whitenessCriteria',{}},@est_checkMVARWhiteness,'Compute residual autocorrelation and variance','cat','Validation Methods','suppress',{'WhitenessCriteria','MultipleComparisonsCorrection','SignificanceLevel','prctWinToSample','VerbosityLevel'}), ...
    arg_subtoggle({'checkConsistency','CheckConsistency'},{'verb',verb,'prctWinToSample',prctWinToSample},@est_checkMVARConsistency,sprintf('Check model consistency. \nA "consistent" model is able to generate data with similar covariance structure as the original data'),'cat','Validation Methods','suppress',{'prctWinToSample','VerbosityLevel'}), ...
    arg_subtoggle({'checkStability','CheckStability'},{'verb',verb},@est_checkMVARStability,sprintf('Check model stability. \nA stable model is one that cannot "blow up" (e.g. generate data that increases to infinity) -- this is a requirement for causal modeling.'),'cat','Validation Methods','suppress',{'prctWinToSample','VerbosityLevel'}), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Data Reduction'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model. Default is empty (use all windows)','cat','Data Reduction'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg({'plot','PlotResults'},true,[],'Plot results') ...
    );

% commit EEG variable to workspace
% [data g] = hlp_splitstruct(g,{'g.EEG'});
% arg_toworkspace(data);
% clear data;

% initialize default output
[whitestats PCstats stabilitystats residualstats] = deal([]);

if ~isfield(g.EEG.CAT,'MODEL')
    error('est_validateMVAR:NoModel','g.EEG.CAT.MODEL must be present in the dataset');
else
    MODEL = g.EEG.CAT.MODEL;
end

% determine which windows to use
if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round(MODEL.winStartTimes*g.EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    g.prctWinToSample = 100;
end

if g.verb==2
    % create waitbar
    waitbarTitle = sprintf('Validating Model %s...', ...
        fastif(isempty(g.EEG.condition),'',['for ' g.EEG.condition]));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle, ...
        'Color', hlp_getNextUniqueColor('reset'), ...
        'CanCancel','on', ...
        'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
    %     wb_cleanup = onCleanup(@() multiWaitbar(waitbarTitle,'Close'));
end

% determine number of checks we'll do (for progbar)
numchecks = sum([g.checkConsistency.arg_selection  ...
    g.checkWhiteness.arg_selection    ...
    g.checkStability.arg_selection    ...
    g.checkResidualVariance.arg_selection]);
curcheck  = 0;

% check residual whiteness
% -------------------------------------------------------------------------
if g.checkWhiteness.arg_selection
    [whitestats residualstats] = est_checkMVARWhiteness(...
                                    'EEG',g.EEG,                        ...
                                    'MODEL',MODEL,                      ...
                                     g.checkWhiteness,                  ...
                                    'winStartIdx',g.winStartIdx,        ...
                                    'prctWinToSample',g.prctWinToSample,...
                                    'verb',g.verb);
    if isempty(whitestats)
        % whiteness checks cancelled
        g.checkWhiteness.arg_selection = false;
        numchecks = numchecks - 1;
        curcheck  = curcheck  - 1;
    else
        winTimes = whitestats.winStartTimes;
    end
    if g.verb==2 && ~updateWaitbar(), return; end
end

% check residual variance
% -------------------------------------------------------------------------
if g.checkResidualVariance.arg_selection
    if ~g.checkWhiteness.arg_selection
        [~, residualstats] = est_checkMVARWhiteness(...
                                'EEG',g.EEG,                            ...
                                'MODEL',MODEL,                          ...
                                 g.checkResidualVariance,               ...
                                'whitenessCriteria',{},                 ...
                                'winStartIdx',g.winStartIdx,            ...
                                'prctWinToSample',g.prctWinToSample,    ...
                                'verb',g.verb);
    end
    % bookeeping
    if isempty(residualstats)
        % whiteness checks cancelled
        g.checkResidualVariance.arg_selection = false;
        numchecks = numchecks - 1;
        curcheck  = curcheck  - 1;
    else
        winTimes = residualstats.winStartTimes;
    end
    if g.verb==2 && ~updateWaitbar(), return; end
end

% check model consistency
% -------------------------------------------------------------------------
if g.checkConsistency.arg_selection
    PCstats = est_checkMVARConsistency(...
                    'EEG',g.EEG,                        ...
                    'MODEL',MODEL,                      ...
                     g.checkConsistency,                ...
                    'winStartIdx',g.winStartIdx,        ...
                    'prctWinToSample',g.prctWinToSample,...
                    'verb',g.verb);
    % bookeeping
    if isempty(PCstats)
        % consistency checks cancelled
        g.checkConsistency.arg_selection = false;
        numchecks = numchecks - 1;
        curcheck  = curcheck  - 1;
    else
        winTimes = PCstats.winStartTimes;
    end
    if g.verb==2 && ~updateWaitbar(), return; end
end

% check model stability
% -------------------------------------------------------------------------
if g.checkStability.arg_selection
    stabilitystats = est_checkMVARStability(...
                        'EEG',g.EEG,                        ...
                        'MODEL',MODEL,                      ...
                         g.checkStability,                  ...
                        'winStartIdx',g.winStartIdx,        ...
                        'prctWinToSample',g.prctWinToSample,...
                        'verb',g.verb);
    % bookeeping
    if isempty(stabilitystats)
        % consistency checks cancelled
        g.checkStability.arg_selection = false;
        numchecks = numchecks - 1;
        curcheck  = curcheck  - 1;
    else
        winTimes = stabilitystats.winStartTimes;
    end
    if g.verb==2 && ~updateWaitbar(), return; end
end


% plot results
% -------------------------------------------------------------------------
if g.plot && numchecks > 0
    if g.checkWhiteness.arg_selection
        whitenessCriteria = g.checkWhiteness.whitenessCriteria;
    else
        whitenessCriteria = {};
    end
    vis_plotModelValidation({whitestats},{PCstats},{stabilitystats}, ...
        'whitenessCriteria',whitenessCriteria,                       ...
        'checkWhiteness',g.checkWhiteness.arg_selection,             ...
        'checkConsistency',g.checkConsistency.arg_selection,         ...
        'checkStability',g.checkStability.arg_selection,             ...
        'conditions',{g.EEG.condition},                              ...
        'windowTimes',winTimes);
end

% % clean up
if g.verb==2
    multiWaitbar(waitbarTitle,'Close');
end


    % subfunction updateWaitbar
    % ---------------------------------------------------------------------
    function nocancel = updateWaitbar()
        nocancel = true;
        curcheck = curcheck + 1;
        drawnow;
        cancel = multiWaitbar(waitbarTitle,curcheck/numchecks);
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle, ...
                     'Are you sure you want to cancel all validation?')
            nocancel = false;
            return;
        end
    end

end
