
function [ALLEEG cfg] = pop_est_mvarConnectivity(ALLEEG,typeproc,varargin)
%
% Computes connectivity estimates from a precomputed MODEL
% and stores the result in ALLEEG.CAT.Conn
%
% Input:
%
%   ALLEEG          EEG structure with fields CAT.MODEL
%
% Optional:
%
%   'connmethods':  cell array of strings denoting connectivity methods to
%                   compute (parenthesized acronym from list below)
%
%                     DIRECTED TRANSFER FUNCTION MEASURES:
%                         Directed Tranfer Function (DTF)
%                         Normalized DTF (nDTF)
%                         Direct DTF (dDTF)
%                         Direct DTF (with full causal normalization) (dDTF08)
%                         Full-frequency DTF (ffDTF)
%                      PARTIAL DIRECTED COHERENCE MEASURES
%                         Partial Directed Coherence (PDC)
%                         Normalized PDC (nPDC)
%                         Generalized Partial Directed Coherence (GPDC)
%                         Partial Directed Coherence Factor (PDCF)
%                         Renormalized Partial Directed Coherence (RPDC)
%                      GRANGER-GEWEKE CAUSALITY MEASURES
%                         Granger-Geweke Causality (GGC)
%                      SPECTRAL COHERENCE MEASURES
%                         Complex Coherence (Coh)
%                         Imaginary Coherence (iCoh)
%                         Partial Coherence (pCoh)
%                         Multiple Coherence (mCoh)
%                      SPECTRAL DENSITY MEASURES
%                         Complex Spectral Density (S)
%   'absvalsq':     Boolean (true,false) determining whether to return the
%                   square of the absolute value of complex measures (def:
%                   true)
%   'spectraldecibels': Boolean (true,false) determining whether to return
%                       the spectral power in units of decibels
%                       (10*log10(Power)) (def: false)
%
% Output:
%
%   ALLEEG          EEG structure with results in EEG.CAT.Conn.
%                   Conn.(connmethod) is a [num_chans x num_chans x
%                   num_freqs x num_time] connectivity matrix
%   params          The options used in the connectivity estimation
%
% See Also: est_mvarConnectivity(), pop_est_fitMVAR()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
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
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG(1),'MODEL',ALLEEG(1).CAT.MODEL},varargin]),false);
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

% initialize progress bar
if cfg.verb==2 && length(ALLEEG)>1
    waitbarTitle = 'Estimating Connectivity';
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle,...
                 'Color', [0.8 0.0 0.1],  ...
                 'CanCancel','on',        ...
                 'CancelFcn',@(a,b)disp('[Cancel requested. Please wait...]'));
end

% now calculate connectivity
for cnd=1:length(ALLEEG)
    [Conn] = feval(fcnHandle,'ALLEEG',ALLEEG(cnd),'MODEL',ALLEEG(cnd).CAT.MODEL,cfg);
    if ~isempty(Conn)
        ALLEEG(cnd).CAT.Conn = Conn; 
    end
    % clear any existing visualization GUI config files
    visFields = fieldnames(ALLEEG(cnd).CAT.configs);
    visFields = visFields(~cellfun(@isempty,strfind(visFields,'vis_')));
    for k=1:length(visFields)
        ALLEEG(cnd).CAT.configs.(visFields{k}) = struct([]);
    end
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
    
    if cfg.verb==2 && length(ALLEEG)>1
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,cnd/length(ALLEEG));
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            break;
        end
    end
    
end

% cleanup progress bar
if cfg.verb==2 && length(ALLEEG)>1
    multiWaitbar(waitbarTitle,'Close');
end
