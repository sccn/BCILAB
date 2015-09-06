function ci = stat_computeCI(varargin)
% compute confidence intervals across last dimension of matrix PConn
%
% Inputs:
% 
%       PConn:     connectivity distribution object returned from stat_bootstrap()
%                  or a matrix of distribution values for a single
%                  estimator
%       ci_alpha:     ci_alpha-significance threshold (e.g. ci_alpha=0.05)
%       ci_tail:      'upper': compute only upper ci (lower is set to mean of
%                           distribution)
%                  'lower': compute only lower ci (upper is set to mean of
%                           distribution)
%                  'both': compute upper and lower confidence intervals
%                          (ci_alpha is divided by 2)
% Outputs:
%
%       ci:     the mean of the distribution
%
% See Also: stat_bootstrap()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.8
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

arg_define([0 Inf],varargin, ...
    arg_norep({'PConn','SurrogateDistrib'},mandatory,[],'Surrogate distribution.'), ...
    arg({'ci_alpha','Alpha'},0.05,[0 1],'Significance level. For example, a value of alpha=0.05 will produce (1-alpha)*100 = 95% confidence intervals.'), ...
    arg({'ci_tail','Tail'},'both',{'right','left','both','one'},'Tail. One-tailed (right-tailed) or two-tailed test. ''Right'' gives upper confidence interval. ''Left'' gives lower confidence interval. ''One'' defaults to ''right.'''), ...
    arg({'testMethod','TestMethod'},'percentile',{'percentile','corrected percentile','normal'},sprintf('Comparison method.\nQuantile: Determines the quantile of the difference (A-B) distribution at which zero lies. From this, one derives a p-value for rejecting the hypothesis that paired samples from both conditions are equal.')), ...
    arg_norep({'Conn','Obs'},[],[],'Estimator') ...
    );


if isstruct(PConn)
    % multiple connectivity methods
    connmethods = hlp_getConnMethodNames(PConn);
    for m=1:length(connmethods)
        % compute stats for all methods
        ci.(connmethods{m}) = stat_computeCI(PConn(cnd).(connmethods{m}),ci_alpha,ci_tail);  
    end
else
    % PConn is a matrix 
    sz = size(PConn);
    nd = length(sz);
    ci_alpha = 100*ci_alpha;
    
    switch lower(ci_tail)
        case 'both'
            lb = prctile(PConn,ci_alpha/2,nd);      % lb
            ub = prctile(PConn,100-ci_alpha/2,nd);  % ub
        case {'ub' 'left'}
            mval = mean(PConn,nd); % mean of estimator
            lb = mval;
            ub = prctile(PConn,100-ci_alpha,nd);    % ub
        case {'lb' 'right'}
            mval = mean(PConn,nd); % mean of estimator
            lb = prctile(PConn,ci_alpha,nd);        % lb
            ub = mval;
        otherwise
            error('SIFT:stat_computeCI','unknown ci_tail option');
    end
    
    nd = ndims(lb);
    ci = cat(nd+1,lb,ub);         % append ci to last dim
    ci = permute(ci,[nd+1 1:nd]); % permute to first dim
end

