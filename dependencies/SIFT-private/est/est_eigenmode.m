
function EIGMODE = est_eigenmode(varargin)
%
% Estimate eigenmodes of VAR process. EEG.CAT.MODEL must be present as 
% output by est_fitMVAR with arfit option selected. Wrapper function for 
% armode() from the ARfit package by T. Schneider and A. Neumaier [1]
%
% This requires the function armode() from ARFIT package (Schneider, 2000)
%
% Inputs:
%
%       EEG:        EEG data structure with EEG.CAT.MODEL present
%       verb:       verbosity (true/false)
%
% Outputs:
%
%   EIGMODE
%       .modes:         Columns contain the estimated eigenmodes of the VAR model.  
%       .modeconf:      Margins of error for the components of the estimated 
%                       eigenmodes, such that (S +/- Serr) are approximate 95% 
%                       confidence intervals for the individual components of 
%                       the eigenmodes.
%       .period:        The first row contains the estimated oscillation period
%                       for each eigenmode. The second row contains margins of
%                       error such that ( period(1,k) +/- period(2,k) ) are
%                       approximate 95% confidence intervals for the period of
%                       eigenmode .modes(:,k).  For a purely relaxatory eigenmode, 
%                       the period is infinite (Inf). For an oscillatory eigenmode, 
%                       the periods are finite.
%       .dampingTime:   The first row contains the estimated oscillation
%                       damping time for each eigenmode. The second row contains 
%                       margins of error such that ( dampingTime(1,k) +/- dampingTime(2,k) )
%                       are approximate 95% confidence intervals for the damping time of
%                       eigenmode .modes(:,k)
%       .exctn:         The excitation of an eigenmode measures its dynamical importance
%                       and is returned as a fraction exctn that is normalized such that
%                       the sum of the excitations of all eigenmodes equals one.
%       .lambda:        The columns contain the eigenvalues of the
%                       eigenmodes
%
%
% See Also: est_fitMVAR(), armode()
%
% References:
%
% [1] Schneider T, Neumaier A (2001) Algorithm 808: ARfit---a matlab package
% for the estimation of parameters and eigenmodes of multivariate 
% autoregressive models. ACM Transactions on Mathematical Software 27:58-65
% http://www.gps.caltech.edu/~tapio/arfit/
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Wrapper function for armode() from the ARfit package by T. Schneider and 
% A. Neumaier [1]
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

arg_define([0 2],varargin,...
        arg_norep('EEG',[],[],'EEG data structure. Must contain .CAT.MODEL structure computed using SIFT'), ...
        arg({'verb','Verbosity'},true,[],'Verbosity Level'), ...
        arg({'computeStats','ComputeStats'},true,[],'Compute statistics. This can significatly increase runtime') ...
        );

if ~exist('armode','file')
    error('Unable to locate armode() function. Is ARfit downloaded and in the path?');
end

% starting point of each window (points)
winStartIdx  = floor(EEG.CAT.MODEL.winStartTimes*EEG.srate)+1;    

numWins = length(winStartIdx);

if verb, h=waitbar(0,'performing eigendecomposition ...'); end
if ~computeStats || ~isfield(EEG.CAT.MODEL,'th') || isempty(EEG.CAT.MODEL.th)
    EEG.CAT.MODEL.th = cell(1,numWins); 
end
m = EEG.CAT.nbchan;
for t=1:numWins
    [EIGMODE.modes{t}, EIGMODE.modeconf{t}, EIGMODE.period{t}, EIGMODE.dampingTime{t}, EIGMODE.exctn{t}, EIGMODE.lambda{t}] = armode2(EEG.CAT.MODEL.AR{t}, EEG.CAT.MODEL.PE{t}(:,1:m), EEG.CAT.MODEL.th{t});
    
    if verb, waitbar(t/numWins,h,sprintf('performing eigendecomposition (%d/%d)...',t,numWins)); end
end

if verb, close(h); end