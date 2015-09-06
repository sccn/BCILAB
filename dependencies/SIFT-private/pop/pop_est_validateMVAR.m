function [ALLEEG cfg] = pop_est_validateMVAR(ALLEEG,typeproc,varargin)
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
%     cfg:                    parameter configuration object
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

if nargin<2
    typeproc = 0;
end

fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

% check the dataset
res = hlp_checkeegset(ALLEEG,{'model'});
if ~isempty(res)
    error(['SIFT:' fcnName],res{1});
end

if isfield(ALLEEG(1).CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG(1)},varargin]),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],ALLEEG(1),varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
        cfg = [];
        return;
    end
    
    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

drawnow;

if strcmpi(typeproc,'cfg_only')
    return;
end

% Apply model validation routines
for cnd=1:length(ALLEEG)
    
    [ALLEEG(cnd).CAT.VALIDATION.whitestats ...
     ALLEEG(cnd).CAT.VALIDATION.PCstats    ...
     ALLEEG(cnd).CAT.VALIDATION.stabilitystats ...
     ALLEEG(cnd).CAT.VALIDATION.residualstats] ...
        = est_validateMVAR('EEG',ALLEEG(cnd),cfg);
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
end