function [ALLEEG cfg] = pop_est_selModelOrder(ALLEEG,typeproc,varargin)
%
% Fit a series of MVAR models up to a specified model order and compute the
% model order selection (information) criteria. If the model fitting
% algorithm is 'vierra-morf' then we fit a single model up to maximum order
% and downdate the prediction errors. This function generates a figure
% containing the results of the model fitting.
% If nargin == 2, an input GIU will be generated. Otherwise
% est_selModelOrder is called directly.
%
% Input:
%
%   EEG                Preprocessed EEG structure. Must contain .CAT
%   typeproc           Reserved for future use. Use 0.
%
% Optional:            arguments to est_selModelOrder()
%
%     'icselector'         cell array of strings denoting which model order
%                          selection criteria to estimate
%     'downdate'           [true, false] whether or not to use downdated noise
%                          covariance matrices
%                          ('aic','sbc','fpe','hq')
%     'algorithm',         string denoting which algorithm to use for model
%                           fitting ('vierra-morf','arfit')
%     'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
%     'morder',            [min max] VAR model order to fit
%     'winlen',            window length (sec)
%     'winstep',           window step size (sec)
%     'epochTimeLims',     time range to analyze (sec) where 0 = start of the epoch
%     'prctWinToSample',   percent of time windows to randomly select  [0 100]
%     'verb',              verbosity level (0=no output, 1=text, 2=gui)
%     'normalize'          cell array containing one or more of
%                           {'temporal', 'ensemble'}. This performs ensemble
%                           normalization or temporal normalization (or both)
%                           within each window
%
% Output:
%
%   IC                 a structure containing results of model order selection
%                      IC.selector     - the chosen information criteria
%                      IC.pmin         - the minimum model order tested
%                      IC.pmax         - the maximum model order tested
%                      IC.('sel') contains results for a selector 'sel'.
%                      This consists of subfields
%                           .ic         - [P numwins] matrix of information
%                                         critera for all P model orders tested
%                                         P = morder(2)-morder(1)+1 is the
%                                         number of model orders tested
%                           .minic      - the minimum of ic across model
%                                         orders
%                           .popt       - the model order that minimizes ic
%
% See Also: est_selModelOrder()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%   Springer.
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
res = hlp_checkeegset(ALLEEG,{'cat'});
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

% Apply model selection
for cnd=1:length(ALLEEG)
    % calculate the information criteria
    ALLEEG(cnd).CAT.IC = est_selModelOrder('EEG',ALLEEG(cnd),cfg);
    
    if isempty(ALLEEG(cnd).CAT.IC)
        % use canceled
        return;
    end
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
end

% determine if we will proceed to model fitting using these options
res=questdlg2(sprintf(['Do you want to proceed to model fitting?\n' ...
                       'A Model-fitting GUI will be generated for you based on the options you selected above']), ...
              'Model Order Selection Assistant', 'No', 'Yes', 'Yes');
          
if strcmpi(res,'yes')
    % Open the model-fitting GUI for model fitting. 
    % Once model is fit results will be return in EEG structure
    modFuncName = ['pop_' ALLEEG(1).CAT.IC.modelFitting.modelingFuncName];
    ALLEEG = feval(modFuncName, ALLEEG,0,ALLEEG(1).CAT.IC.modelFitting.modelingArguments);
end



