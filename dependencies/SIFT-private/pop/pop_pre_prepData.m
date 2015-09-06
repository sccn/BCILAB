function [ALLEEG_out cfg] = pop_pre_prepData(ALLEEG,typeproc,varargin)
%
% Preprocess EEG dataset(s) for connectivity analysis. See [1] for
% mathematical details on preprocessing steps.
%
%
% Input:
%
%   ALLEEG:         Array of EEGLAB datasets to preprocess.
%   typeproc:       Reserved for future use. Use 0
%
% Optional:
%
%   <'Name',value> pairs as defined in pre_prepData()
%
% Output:
%
%   ALLEEG:         Prepocessed EEG structure(s)
%   cfg:            Argument specification structure.
%
%
% See Also: pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2009, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu
%
% Revised Jan 2010.

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

% set default output
ALLEEG_out = ALLEEG;
cfg = [];
        
% generate splash screen
% initialize SIFT, etc
StartSIFT(~strcmpi(typeproc,'nogui'),true);

fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

% check if we've applied SIFT to this dataset before
res = hlp_checkeegset(ALLEEG,{'cat'});
if isempty(res) && isfield(ALLEEG(1).CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the default config from function and overload supplied args
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG(1)},varargin]),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],ALLEEG(1),varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
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

% initialize progress bar
if cfg.verb==2 && length(ALLEEG)>1
    waitbarTitle = 'Preprocessing datasets';
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle,...
                 'Color', [0.8 0.0 0.1],  ...
                 'CanCancel','on',        ...
                 'CancelFcn',@(a,b)disp('[Cancel requested. Please wait...]'));
end


% re-initialize output
clear ALLEEG_out;

% preprocess datasets
for cnd=1:length(ALLEEG)
    
    % execute the low-level function
    [ALLEEG_out(cnd)] = feval(fcnHandle,'EEG',ALLEEG(cnd),cfg);
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG_out(cnd).CAT.configs.(fcnName) = cfg;
    end

    
    if cfg.verb==2 && length(ALLEEG)>1
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,cnd/length(ALLEEG));
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            % restore original dataset
            ALLEEG_out = ALLEEG;
            break;
        end
    end
    
        
end

% cleanup progress bar
if cfg.verb==2 && length(ALLEEG)>1
    multiWaitbar(waitbarTitle,'Close');
end

