function [EEGout cfg] = pop_sim_eegdata(EEG,typeproc,varargin)
%
% Simulate EEG data using a dynamical modeling framework and 
% a forward head model.
%
% Input:
% Optional:
%
%   EEG:            existing EEG dataset (configs will be retrived from
%                   here)
%   typeproc:       if 'nogui' don't generate GUI
%
%   <'Name',value> pairs as defined in sim_varmodel()
%
% Output:
%
%   EEGout:         Simulated EEG structure(s).
%                   Optionally this may be an array of
%                   EEG structs with the second struct 
%                   being the ground truth model.
%   cfg:            Argument specification structure.
%
%
% See Also: sim_varmodel(), sim_eegdata(),
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu
%

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

if nargin<1
    EEG = [];     
end
if nargin<2
    typeproc = 0; 
end
EEGout = EEG;
cfg    = [];
        
fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

if isempty(hlp_checkeegset(EEG,{'cat'})) && isfield(EEG.CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(EEG.CAT.configs.(fcnName)) varargin];
end

if strcmpi(typeproc,'nogui')
    % get the default config from function and overload supplied args
    cfg = arg_tovals(arg_report('rich',fcnHandle,varargin),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
        return;
    end
    
    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

drawnow;

if ~cfg.srcdyn.makeEEGset.arg_selection
    error('SIFT:sim_varmodel',['If using pop_' fcnName '(), you must enable the BuildEEGLABStructure option.\n' ...
                               'Use ' fcnName '() from the command-line to return a raw dataset']);
end

if strcmpi(typeproc,'cfg_only')
    return;
end

% execute the low-level function
[EEGout GroundTruth] = feval(fcnHandle,cfg);

if ~isempty(GroundTruth)
    EEGout = eeg_store(EEGout,GroundTruth,2);
end
if ~isempty(cfg)
    for k=1:length(EEGout)
        % store the configuration structure
        EEGout(k).CAT.configs.(fcnName) = cfg;
    end
end


