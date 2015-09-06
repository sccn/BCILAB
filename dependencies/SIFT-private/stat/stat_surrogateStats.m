function [Stats ConnNew cfg] = stat_surrogateStats(varargin)
%
% Return surrogate statistics based on a surrogate distribution (bootstrap, jacknife, etc).
% Each surrogate distribution should approximate the distribution of the
% estimator.
%
%
% Stats are performed across the 5th dimension (boostrap samples or subjects).
% This function calls the statcond() routine written by Arnaud Delorme.
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for
%    updated version
% ===============================================
%
% Input                         Information
% ----------------------------------------------------------------------------------------------------------------
% BootstrapConnectivity         Surrogate connectivity structure as returned
%                               by stat_surrogate()
%
% Optional                      Information
% ----------------------------------------------------------------------------------------------------------------
% NullDistribution              Null distribution structure as returned by
%                               stat_surrogate()
%
%
% ConnectivityMethods:          Connectivity estimator(s) to bootstrap
%                               All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn
%                               Possible values: ''
%                               Default value  : 'n/a'
%                               Input Data Type: boolean
%
% MultipleComparisonCorrection: Correction for multiple comparisons.
%                               'numvars' does a bonferonni correction
%                               considering M^2 indep. degrees of freedom,
%                               where M is the dimension of the VAR model
%                               (i.e. number of channels)
%                               Possible values: 'none','fdr','bonferonni','numvars'
%                               Default value  : 'fdr'
%                               Input Data Type: string
%
% ComputeConfIntervals:         Compute confidence intervals
%                               Input Range  : Unrestricted
%                               Default value: 0
%                               Input Data Type: boolean
%
%     | Alpha:                  Significance level (two sided)
%                               A value of 0.95 will produce 95% confidence intervals.
%                               Input Range  : [0  1]
%                               Default value: 0.95
%                               Input Data Type: real number (double)
%
% statcondargs:                 List of paired arguments for statcond()
%                               Possible values: Unrestricted
%                               Default value  : 'mode','perm'
%                               Input Data Type: any evaluable Matlab expression.
%
% VerbosityLevel:               Verbosity level. 0 = no output, 1 = verbose output
%                               Possible values: 0,1
%                               Default value  : 1
%                               Input Data Type: real number (double)
%
% Output                        Information
% ----------------------------------------------------------------------------------------------------------------
% Stats                         Contains p-values, confidence interval, and
%                               other statistics for each connectivity
%                               estimator in BootstrapConnectivity. The
%                               dimensions for the stats matrices are generally
%                               the same (or +1) as the dimensions of the matrix
%                               for the original estimator.
%                               Output format is as follows:
%                               Stats.<estimator>.[pval, ci, thresh, ...]
% ConnMean                      Contains the mean of the bootstrap
%                               distribution for each connectivity
%                               estimator (if using bootstrap methods).
%                               This does not apply to PhaseRand methods
%
%
% See Also: stat_surrogate(), stat_analyticStats(), statcond()
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 5.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD.
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



% extract some stuff from inputs for arg defaults
EEG = arg_extract(varargin,{'EEG','ALLEEG'},1);
if ~isempty(EEG(1)) && isempty(hlp_checkeegset(EEG,{'pconn'}))
    for k=1:length(EEG)
        PConn(k) = EEG(k).CAT.PConn;
    end
    ConnNames   = hlp_getConnMethodNames(PConn(1));
    conndef     = ConnNames;
else
    ConnNames = {''};
    conndef = '';
    PConn   = [];
end


% statTestDef = {'Hnull','Hbase','Hab'};

g = arg_define([0 Inf],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEG structure. This must contain a surrogate distribution in EEG.CAT.PConn'), ...
    arg_subswitch({'statTest','StatisticalTest'},hlp_getAvailableStatTests(EEG,true),...
    hlp_getAvailableStatTests(EEG,false),{'Statistical test to perform.', ...
       sprintf(['\n' ...
                'Hnull: Compare PConn to null distribution (must be provided in ''NullDistribution'' argument) (one-tailed, unpaired).\n' ...
                'Hbase: Compare each sample in PConn to baseline (one- or two-tailed, unpaired).\n' ...
                'Hab:   Compute significant difference for two conditions (PConn(1)-PConn(2)) (one- or two-tailed, paired).'])}), ...  % ,'suppress','Hbase'
    arg({'connmethods','ConnectivityMethods'},conndef,ConnNames,'Connectivity estimator(s) to bootstrap. All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn','type','logical'), ...
    arg({'verb','VerbosityLevel'},1,{int32(0) int32(1)},'Verbosity level. 0 = no output, 1 = verbose output') ...
    );

[Stats ConnNew] = deal([]);

% check if we have surrogate distributions
% ... or null distribution
res = {hlp_checkeegset(g.EEG,{'pconn'})
       hlp_checkeegset(g.EEG,{'pnull'})};
idx = ~cellfun(@isempty,res);
if all(idx==1)
    error(res{idx}{1});
end
if ~ismember_bc('mode',lower(g.statTest.statcondargs))
    g.statTest.statcondargs = [g.statTest.statcondargs 'mode', 'perm'];
end
if g.verb==0
    g.statTest.statcondargs = [g.statTest.statcondargs 'verbose', 'off'];
else
    g.statTest.statcondargs = [g.statTest.statcondargs 'verbose', 'on'];
end
if nargout > 2, cfg = g; end

if ~exist('PConn','var')
    for k=1:length(EEG)
        PConn(k) = EEG(k).CAT.PConn; 
    end
end

for m=1:length(g.connmethods)
    
    switch g.statTest.arg_selection
        case 'Hnull'
            % Our null hypothesis is Conn(i,j)=0
            % perform one-sided statistical test against null hypothesis
            
            Stats.tail = g.statTest.tail;
            
            if length(PConn)>1
                error('SIFT:stat_surrogateStats','Please select a single condition for Hnull test.');
            end
            
            if strcmpi(PConn.mode,'phaserand')
                % We are testing with respect to a phase-randomized null
                % distribution. A p-value for rejection of the null hypothesis
                % can be obtained by computing the probability that the
                % observed connectivity is a random sample from the null distribution
                
                if g.verb
                    fprintf('\nTesting estimator %s\n',g.connmethods{m});
                    fprintf('Computing statistics against null hypothesis C(i,j)=0\n');
                    fprintf('Stats are based on %s distribution\n',PConn.mode);
                end
                
%                 if ndims(PConn.(g.connmethods{m}))==ndims(g.null.(g.connmethods{m}))
%                     % a bootstrap distribution for the estimator was provided so
%                     % we need to compute the mean of the bootstrap
%                     % distribution for comparison to null distribution
%                     PConn.(g.connmethods{m}) = stat_getDistribMean(PConn.(g.connmethods{m}));
%                 end
                
                switch g.statTest.testMethod
                    case 'quantile'
                        Stats.tail = 'right';   
                        % get empirical p-values against the null hypothesis that
                        % observed data comes from the null distribution
                        Stats.(g.connmethods{m}).pval = stat_surrogate_pvals( ...
                                        PConn.(g.connmethods{m}),             ...
                                        g.EEG.CAT.Conn.(g.connmethods{m}),    ...
                                        Stats.tail);
                        % compute the threshold based on 1-alpha percentile of null
                        % distribution
                        Stats.(g.connmethods{m}).thresh = prctile( ...
                                        PConn.(g.connmethods{m}),  ...
                                        100-100*g.statTest.alpha,  ...
                                        ndims(PConn.(g.connmethods{m})));
                
                        % there are no confidence intervals to estimate here
                        Stats.(g.connmethods{m}).ci = [];
                        
                    case 'statcond'
                        [statval, df, Stats.(g.connmethods{m}).pval] = statcond( ...
                            {PConn.(g.connmethods{m})},'mode','perm', ...
                            'tail',Stats.tail,g.statTest.statcondargs{:});
                        Stats.(g.connmethods{m}).ci = [];
                end
            else
                % we are testing w.r.t. the approximate distribution of
                % the estimator itself. A p-value for rejection of the null
                % hypothesis can be obtained by computing the probability
                % that a sample from the estimator's distribution is less
                % than or equal to zero
                % univariate nonparametric test
                if g.verb
                    fprintf('\nTesting estimator %s\n',g.connmethods{m});
                    fprintf('Computing statistics against null hypothesis C(i,j)=0\n');
                    fprintf('Stats are based on %s distribution\n',PConn.mode);
                end

                switch g.statTest.testMethod
                    case 'quantile'
                    	error('You must first compute a null (e.g. phase-randomized) distribution. See pop_stat_surrogateGen()');
                    case 'statcond'
                        [statval, df, Stats.(g.connmethods{m}).pval] = statcond( ...
                            {PConn.(g.connmethods{m})},'mode','perm', ...
                            'tail',Stats.tail,g.statTest.statcondargs{:});
                end
            end

        case 'Hab'
            % For conditions A and B, the null hypothesis is either
            % A(i,j)<=B(i,j), for a one-sided test, or
            % A(i,j)=B(i,j), for a two-sided test
            % A p-value for rejection of the null hypothesis can be
            % obtained by taking the difference of the distributions
            % and computing the probability
            % that a sample from the difference distribution is non-zero
            % This is a paired nonparametric test
            
            Stats.tail = g.statTest.tail;
            
            if length(PConn)~=2
                error('SIFT:stat_surrogateStats','BootstrapConnectivity object must have two conditions for difference test');
            elseif size(PConn(1).(g.connmethods{m}))~=size(PConn(2).(g.connmethods{m}))
                error('SIFT:stat_surrogateStats','BootstrapConnectivity matrices must be of equal dimension for both conditions');
            end
            
            if g.verb
                fprintf('\nTesting estimator %s\n',g.connmethods{m});
                fprintf('Computing statistics against null hypothesis A(i,j)=B(i,j)\n');
                fprintf('This is a %s-sided test for significant differences between conditions\n',Stats.tail);
                fprintf('Stats are based on %s distributions\n',PConn(1).mode);
            end
            
            % get the order in which to subtract the datasets...
            [dummy setIdx] = ismember_bc(g.statTest.datasetOrder,hlp_getCondOrderings(EEG));
            Stats.diffOrder = fastif(setIdx==1,[1 2],[2 1]);
            % ... and reorder datasets
            PConn = PConn(Stats.diffOrder);
            
            switch g.statTest.testMethod
                case 'quantile'
                    sz = size(PConn(1).(g.connmethods{m}));
                    Pdiff = PConn(1).(g.connmethods{m}) ...
                          - PConn(2).(g.connmethods{m});
                    Stats.(g.connmethods{m}).pval = stat_surrogate_pvals(...
                                        Pdiff,zeros(sz(1:end-1)),Stats.tail);
                    %[statval, df, Stats.(g.connmethods{m}).pval] = statcond( { }, 'mode','perm', ..
                    % 'surrog', Pdiff, 'stats', zeros(sz(1:end-1)), ..
                    % 'tail',Stats.tail,g.statTest.statcondargs{:});
%                     Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval;

                    if g.statTest.computeci
                        % compute two-sided confidence intervals
                        Stats.(g.connmethods{m}).ci = stat_computeCI( ...
                                        Pdiff,g.statTest.alpha,Stats.tail);
                    end
                    
                    if nargout>1
                        % return the expected value of the difference
                        % between conditions. Result is stored in ConnNew
                        if ~exist('ConnNew','var') || isempty(ConnNew)
                            ConnNew = rmfield(PConn(1), ...
                                      ['resampleTrialIdx',hlp_getConnMethodNames(PConn)]);
                        end
                        ConnNew.(g.connmethods{m}) = mean(Pdiff,ndims(Pdiff));
                    end
                 case 'statcond'
                    [statval, df, Stats.(g.connmethods{m}).pval] = statcond(...
                            {PConn.(g.connmethods{m})},'mode','perm', ...
                            'tail',Stats.tail,g.statTest.statcondargs{:});
            end
            
        case 'Hbase'
            % For conditions A, the null hypothesis is
            % C(i,j)=baseline_mean(C). This is a two-sided test.
            % A p-value for rejection of the null hypothesis can be
            % obtained by obtaining the distribution of the difference from
            % baseline mean and computing the probability
            % that a sample from this distribution is non-zero
            % This is a paired nonparametric test
            
            Stats.tail = g.statTest.tail;
            
            if length(g.EEG)>1
                error('SIFT:stat_surrogateStats','Please select a single condition for Hbase test.');
            end
            
            if length(PConn.winCenterTimes)==1
                error('SIFT:stat_surrogateStats','Estimator must be time-varying to compute deviation from temporal baseline');
            end
            
            if g.verb
                fprintf('\nTesting estimator %s\n',g.connmethods{m});
                fprintf('Computing statistics against null hypothesis C(i,j)=baseline_mean(C)\n');
                fprintf('This is a %s-sided test for significant deviation from baseline\n',Stats.tail);
                fprintf('Stats are based on %s distributions\n',PConn.mode);
            end
            
            switch g.statTest.testMethod
                case 'quantile'
                    % quantile test
                    % compute the baseline deviation distribution
                    % (Conn - baseline_mean(Conn))
                    if ~g.statTest.testMeans
                        error('TestMeans must be enabled when using ''quantile'' method');
                    end
                    Pdiff = PConn.(g.connmethods{m}) ...
                          - stat_getBaselineDistrib(PConn.(g.connmethods{m}), ...
                                 g.statTest.baseline,PConn.erWinCenterTimes);
                    sz = size(Pdiff);
                    Stats.(g.connmethods{m}).pval = stat_surrogate_pvals( ...
                                    Pdiff,zeros(sz(1:end-1)),Stats.tail);
                    % compute two-sided confidence intervals
                    Stats.(g.connmethods{m}).ci = stat_computeCI( ...
                                    Pdiff,g.statTest.alpha,Stats.tail);
                    
                    if nargout>1
                        % return the expected value of the difference
                        % from baseline mean. Result is stored in ConnNew
                        if ~exist('ConnNew','var') || isempty(ConnNew)
                            ConnNew = rmfield(PConn, ...
                                      ['resampleTrialIdx',hlp_getConnMethodNames(PConn)]);
                        end
                        ConnNew.(g.connmethods{m}) = mean(Pdiff,ndims(Pdiff));
                    end
                case 'statcond'
                    Porig = PConn.(g.connmethods{m});
                    Pbase = stat_getBaselineDistrib(...
                                PConn.(g.connmethods{m}), ...
                                g.statTest.baseline,      ...
                                PConn.erWinCenterTimes,   ...
                                ~g.statTest.testMeans);        
                    [statval, df, Stats.(g.connmethods{m}).pval] = statcond(...
                        {Porig Pbase},'mode','perm', ...
                        'tail',Stats.tail,g.statTest.statcondargs{:});
            end
            Stats.baseline = g.statTest.baseline;
            
            
        case 'BasicStats'
            % compute confidence intervals and means from boostrap or
            % jacknife distributions of the estimator
                        
            if length(g.EEG)>1
                error('SIFT:stat_surrogateStats','Please select a single condition for Hbase test.');
            end
            
            if g.verb
                fprintf('\nTesting estimator %s\n',g.connmethods{m});
                fprintf('Computing confidence intervals and means\n');
                fprintf('Stats are based on %s distributions\n',PConn.mode);
            end
           
            
            if g.statTest.computeci.arg_selection
                Stats.tail = g.statTest.computeci.tail;

                switch g.statTest.computeci.testMethod
                
                    case 'quantile'
                        % quantile test
                        % compute the baseline deviation distribution
                        % (Conn - baseline_mean(Conn))
                        
                        % compute two-sided confidence intervals
                        Stats.(g.connmethods{m}).ci = stat_computeCI(   ...
                                            PConn.(g.connmethods{m}),   ...
                                            g.statTest.computeci.alpha, ...
                                            Stats.tail);
                    case 'stderr'
                        %alpha = g.statTest.computeci.alpha;
                        % return confidence intervals based on the standard error
                        %critval = norminv([alpha/2 1-alpha/2],0,1)
                        
                    case 'statcond'
                end
            end
            
            % ugly hack!
%             sz = size(PConn.(g.connmethods{m}));
%             Stats.(g.connmethods{m}).pval = zeros(sz(1:end-1),'single');
            g.statTest.alpha = g.statTest.computeci.alpha;
            
            % compute mean of the estimator
            if g.statTest.computemean && nargout>1
                % return the expected value of the surrogate distribution. 
                % Result is stored in ConnNew
                if ~exist('ConnNew','var') || isempty(ConnNew)
                    ConnNew = rmfield(PConn, ...
                              ['resampleTrialIdx',hlp_getConnMethodNames(PConn)]);
                end
                nd = ndims(PConn.(g.connmethods{m}));
                ConnNew.(g.connmethods{m}) = mean(PConn.(g.connmethods{m}),nd);
            end
            
    end
    
    % Correct for multiple comparisons
    switch g.statTest.mcorrection
        case 'fdr'
            Stats.(g.connmethods{m}).pval = fdr(Stats.(g.connmethods{m}).pval);
        case 'bonferonni'
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval * numel(Stats.(g.connmethods{m}).pval);
        case 'numvars'
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval * size(Stats.(g.connmethods{m}).pval,1)^2;
        case 'custom_factor'
            Stats.(g.connmethods{m}).pval = Stats.(g.connmethods{m}).pval * custom_factor;
    end
    
    % check if a singleton dimension was squeezed out and, if so,
    % restore the singleton dim
    
    % check pval
    if isfield(Stats.(g.connmethods{m}),'pval')
        szp = size(PConn(1).(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).pval);
        [dummy dimidx] = setdiff_bc(szp(1:end-1),szs);
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).pval = hlp_insertSingletonDim(Stats.(g.connmethods{m}).pval,dimidx+1);
        end
    end
    
    % check ci
    if isfield(Stats.(g.connmethods{m}),'ci')
        szp = size(PConn(1).(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).ci);
        [dummy dimidx] = setdiff_bc(szp(1:end-1),szs(2:end));
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).ci = hlp_insertSingletonDim(Stats.(g.connmethods{m}).ci,dimidx+1);
        end
    else
        Stats.(g.connmethods{m}).ci = [];
    end
    
    % check thresh
    if isfield(Stats.(g.connmethods{m}),'thresh')
        szp = size(PConn(1).(g.connmethods{m}));
        szs = size(Stats.(g.connmethods{m}).thresh);
        [dummy dimidx] = setdiff_bc(szp(1:end-1),szs);
        if ~isempty(dimidx)
            % a singleton dimension was squeezed out, restore it
            Stats.(g.connmethods{m}).thresh = hlp_insertSingletonDim(Stats.(g.connmethods{m}).thresh,dimidx+1);
        end
    else
%         Stats.(g.connmethods{m}).thresh = [];
    end
    
end

statcondargs     = hlp_varargin2struct(g.statTest.statcondargs);
Stats.mode       = statcondargs.mode;
Stats.correction = g.statTest.mcorrection;
Stats.alpha      = g.statTest.alpha;


% helper functions for defining allowable runtime arguments
% -------------------------------------------------------------------------

function arglist = hlp_getAvailableStatTests(EEG,defNameOnly)
% return a cell array (arglist) defining the available available statistical
% tests (null hypotheses) for a given surrograte distribution
% each entry of the cell array is a cell array of the form
%
% hlp_getAvailableStatTests(...,true) returns only the name of the 
%  default modeling approach (as a string)

if isempty(EEG)
    arglist = {};
    return;
end

if nargin<2
    defNameOnly = false; end

if length(EEG)>2
    error('SIFT:stat_surrogateStats','A maximum of two datasets can be compared statistically'); end

% determine available null hypotheses
arglist = {};
if length(EEG)==1
    if strcmpi(EEG.CAT.PConn.mode,'phaserand')
        arglist{end+1} = {'Hnull' @tst_Hnull};
    else
        arglist{end+1} = {'Hbase' @tst_Hbase};  
        arglist{end+1} = {'BasicStats' @tst_basicStats};
    end
end

if length(EEG)==2
    CondDiffOrders = hlp_getCondOrderings(EEG);
    arglist{end+1} = {'Hab' @(varargin) tst_Hab(varargin{:},'dummy',CondDiffOrders)};
end

% arglist = {arglist};

if defNameOnly
    arglist = arglist{1}{1};
    return;
end


function CondDiffOrders = hlp_getCondOrderings(EEG)
% return possible orderings for condition difference
% (A-B), (B-A)

% get the condition names
conditionNames = {EEG.condition};
setNames       = {EEG.setname};
fileNames      = {EEG.filename};

% set up condition difference order defaults
if ~any(cellfun(@isempty,conditionNames))
    % use condition names for labeling
    CondDiffOrders = ...
        {sprintf('%s-%s',conditionNames{1},conditionNames{2}), ...
        sprintf('%s-%s',conditionNames{2},conditionNames{1})};
elseif ~any(cellfun(@isempty,setNames))
    % use set names for labeling
    CondDiffOrders = ...
        {sprintf('%s-%s',setNames{1},setNames{2}), ...
        sprintf('%s-%s',setNames{2},setNames{1})};
elseif ~any(cellfun(@isempty,fileNames))
    % use set filenames for labeling
    CondDiffOrders = ...
        {sprintf('%s-%s',fileNames{1},fileNames{2}), ...
        sprintf('%s-%s',fileNames{2},fileNames{1})};
else
    % use set numbers for labeling
    CondDiffOrders = {'Set 1 - Set 2','Set 2 - Set 1'};
end
    

function args = tst_Hnull(varargin)
args = arg_define(0,varargin, ...
        arg({'testMethod','TestMethod'},'quantile',{'quantile'},sprintf('Comparison method.\nQuantile: Determines the quantile in which an observed sample lies within the null distribution. From this, one derives a p-value for rejecting the hypothesis that the data comes from the null distribution.')), ...
        arg({'tail','Tail'},'right',{'right','left','both','one'},'Tail. One-tailed (right-tailed) or two-tailed test. Right tailed test gives probability that A > B'), ...
        arg({'alpha','Alpha'},0.05,[0 1],'Significance level. This is used for significance thresholds. For example, a value of alpha=0.05 will produce p < alpha=0.05 threshold.'), ...
        arg({'mcorrection','MultipleComparisonCorrection'},'fdr',{'none','fdr','bonferonni','numvars','custom_factor'},'Correction for multiple comparisons. Note: ''numvars'' does a bonferonni correction considering M^2 indep. degrees of freedom, where M is the dimension of the VAR model (i.e. number of channels)'), ...
        arg_nogui('statcondargs',{'mode','perm'},{},'List of paired arguments for statcond()','type','expression','shape','row') ...
        );

    
function args = tst_Hbase(varargin)
args = arg_define(0,varargin, ...
        arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point.','shape','row','type','denserealdouble'), ...
        arg({'testMeans','TestMeans'},true,[],sprintf('Test against mean of baseline. If true, the null hypothesis is that an observation is equal to the baseline mean. In this case, we temporally average samples in the baseline window to obtain the (null) distribution of the baseline mean. \nIf false, the null hypothesis is that an observation is equal to any sample in the baseline. In this case, we concatenate all samples in the baseline window to obtain our null distribution.')), ...
        arg({'testMethod','TestMethod'},'quantile',{'quantile'},sprintf('Comparison method.\nQuantile: Determines the quantile of the baseline (mean) distribution at which an observed sample lies. From this, one derives a p-value for rejecting the hypothesis that the sample comes from the baseline (or is equal to the baseline mean).')), ...
        arg({'tail','Tail'},'both',{'right','left','both','one'},'Tail. One-tailed (right-tailed) or two-tailed test. Right tailed test gives probability that A > B'), ...
        arg({'computeci','ConfidenceIntervals'},true, [],'Compute empirical confidence intervals.'), ...
        arg({'alpha','Alpha'},0.05,[0 1],'Confidence interval significance level. For example, a value of alpha=0.05 will produce (1-alpha)*100 = 95% confidence intervals.'), ...
        arg({'mcorrection','MultipleComparisonCorrection'},'fdr',{'none','fdr','bonferonni','numvars','custom_factor'},'Correction for multiple comparisons. Note: ''numvars'' does a bonferonni correction considering M^2 indep. degrees of freedom, where M is the dimension of the VAR model (i.e. number of channels)'), ...
        arg_nogui('statcondargs',{'mode','perm'},{},'List of paired arguments for statcond()','type','expression','shape','row') ...
        );
    

function args = tst_Hab(varargin)
CondDiffOrders = arg_extract(varargin,'dummy',[],[]);
if isempty(CondDiffOrders)
    CondDiffOrders = {CondDiffOrders}; 
end
args = arg_define(0,varargin, ...
        arg({'datasetOrder','DatasetOrder'},CondDiffOrders{1},CondDiffOrders,'Dataset ordering. A - B tests whether A > B.'), ...
        arg({'testMethod','TestMethod'},'quantile',{'quantile'},sprintf('Comparison method.\nQuantile: Determines the quantile of the difference (A-B) distribution at which zero lies. From this, one derives a p-value for rejecting the hypothesis that paired samples from both conditions are equal.')), ...
        arg({'tail','Tail'},'both',{'right','left','both','one'},'Tail. One-tailed (right-tailed) or two-tailed test. Right tailed test gives probability that A > B'), ...
        arg({'computeci','ConfidenceIntervals'},true, [],'Compute empirical confidence intervals.'), ...
        arg({'alpha','Alpha'},0.05,[0 1],'Confidence interval significance level. For example, a value of alpha=0.05 will produce (1-alpha)*100 = 95% confidence intervals.'), ...
        arg({'mcorrection','MultipleComparisonCorrection'},'fdr',{'none','fdr','bonferonni','numvars','custom_factor'},'Correction for multiple comparisons. Note: ''numvars'' does a bonferonni correction considering M^2 indep. degrees of freedom, where M is the dimension of the VAR model (i.e. number of channels)'), ...
        arg_nogui('statcondargs',{'mode','perm'},{},'List of paired arguments for statcond()','type','expression','shape','row'), ...
        arg_norep('dummy',{{}},[],'dummy') ...
        );

    
function args = tst_basicStats(varargin)
args = arg_define(0,varargin, ...
        arg_subtoggle({'computeci','ConfidenceIntervals'},'on', ...
        { ...
        arg({'testMethod','TestMethod'},'quantile',{'quantile'},sprintf('Test method. \n Quantile: compute asymmetric percentile confidence intervals. This is the better choice if your distribution is derived from a bootstrap or if your distribution is not symmetric.\n stderr: use standard error to derive symmetric confidence intervals. This is suitable only for jacknife or leave-k-out cross-validation.'),'cat','ConfidenceInterval'), ...
        arg({'tail','Tail'},'both',{'both'},'Tail. Upper and lower confidence intervals will be returned','cat','ConfidenceInterval'), ...
        arg({'alpha','Alpha'},0.05,[0 1],'Significance level. For example, a value of alpha=0.05 will produce (1-alpha)*100 = 95% confidence intervals.','cat','ConfidenceInterval'), ...
        },'Compute percentile confidence intervals.','cat','ConfidenceInterval'), ...
        arg_nogui({'mcorrection','MultipleComparisonCorrection'},'none',{'none'},'Correction for multiple comparisons. Note: ''numvars'' does a bonferonni correction considering M^2 indep. degrees of freedom, where M is the dimension of the VAR model (i.e. number of channels)'), ...
        arg({'computemean','Mean'},true,[],'Return mean of distribution. If you provided the boostrap, jacknife, or cross-validation distribution, this returns a plug-in for the mean of the estimator. Note the mean is returned in second function output ConnNew'), ...
        arg_nogui('statcondargs',{'mode','perm'},{},'List of paired arguments for statcond()','type','expression','shape','row') ...
        );
