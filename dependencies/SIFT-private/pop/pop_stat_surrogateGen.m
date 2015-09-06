function [ALLEEG cfg] = pop_stat_surrogateGen(ALLEEG,typeproc,varargin)
%
% Return surrogate statistics based on a surrogate distribution (bootstrap, jacknife, etc).
% Each surrogate distribution should approximate the distribution of the
% estimator. 
%
%
% Input:
%
%   ALLEEG:         EEGLAB dataset to preprocess.
%   typeproc:       Reserved for future use. Use 0
%
% Optional:         
%
%   <'Name',value> pairs as defined in stat_surrogate()
%   
% Output:
%
%   ALLEEG:         EEG structure(s) with surrogate distribution object stored in
%                   ALLEEG.CAT.PConn. To obtain statistics, supply this
%                   ALLEEG structure as input to pop_stat_surrogateStats()
%   cfg:            Argument specification structure.
%
%
% See Also: stat_surrogateGen(), pop_stat_surrogateStats(),
%           stat_surrogateStats()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% 
% Author: Tim Mullen 2010-2011, SCCN/INC, UCSD. 
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

if isfield(ALLEEG(1).CAT.configs,fcnName) ...
   && ~isempty(ALLEEG(1).CAT.configs.(fcnName))
    % We have a stat_surrogate() configuration structure from prior use
    % However, the user might have since changed the modeling approach
    % and may instead want to use modeling configs instead
    
    % default case
    useSurrogateConfigDefault = true;
    
    % First, check if the user has fit a model
    res = hlp_checkeegset(ALLEEG,{'model'});
    if isempty(res)
        % the user has fit a model...
        % ...so get the name of the modeling approach function
        modelApproach = ALLEEG(1).CAT.MODEL.modelapproach;
        modFcnName  = hlp_getModelingApproaches('mfileNameOnly',modelApproach);
        
        % ...and check if we have a configuration structure for it
        if ~isempty(ALLEEG(1).CAT.configs.(modFcnName))
            
            % We have a modeling config struct, so check if the user 
            % prefers to ignore the stat_surrogate() configs and instead 
            % set the default options based on the saved configs for 
            % model-fitting, etc.
            res = questdlg2(...
                    sprintf(['It looks like you''ve previously estimated a surrogate distribution.\n', ...
                             'I found a saved set of default options for stat_surrogate()\n', ...
                             'However, it also looks like you have previously fit a model.\n', ...
                             'I can populate the default options using either the stored options for\n', ...
                             'stat_surrogate() or those from modeling configuration.']), ...
                	'Load Default Options', ...
                    'Use Modeling Defaults',...
                    'Use SurrogateStats Defaults', ...
                    'Use Modeling Defaults');
            if strcmpi(res,'Use Modeling Defaults')
                useSurrogateConfigDefault = false;
            end
        end
    end
    
    if useSurrogateConfigDefault
        % get default configuration (from prior use) and merge with varargin
        varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
    end
    
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

% if length(ALLEEG)>1
%     fprintf('More than one dataset is loaded. I will use the same resampling schedule for all datasets\n');
%     fprintf('This will enable us to form the joint bootstrap distribution over multiple conditions\n');
% end

% compute surrogate distributions
for cnd=1:length(ALLEEG)
        
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
   
    % execute the low-level function
    PConn = feval(fcnHandle,'EEG',ALLEEG(cnd),cfg);
    
    if ~isempty(PConn)
        ALLEEG(cnd).CAT.PConn = PConn;
    else
        % user canceled
        return;
    end
    
    % reset the surrogate stats config
    ALLEEG(cnd).CAT.configs.stat_surrogateStats = struct([]);
end