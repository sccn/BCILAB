function [ALLEEG cfg] = pop_stat_surrogateStats(ALLEEG,typeproc,varargin)
%
% Return surrogate statistics based on a surrogate distribution (bootstrap, jacknife, etc).
% Each surrogate distribution should approximate the distribution of the
% estimator. 
%
%
% Input:
%
%   ALLEEG:         EEGLAB dataset to preprocess. Must contain .CAT.PConn
%                   with surrogate distributions
%   typeproc:       Reserved for future use. Use 0
%
% Optional:         
%
%   <'Name',value> pairs as defined in stat_surrogate()
%   
% Output:
%
%   ALLEEG:         EEG structure(s) with Stats object stored in
%                   ALLEEG.CAT.Stats
%   cfg:            Argument specification structure.
%
%
% See Also: stat_surrogate()
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
res = hlp_checkeegset(ALLEEG,{'pconn'});
if ~isempty(res)
    error(['SIFT:' fcnName],res{1});
end

% if isfield(ALLEEG(1).CAT.configs,fcnName)
%     % get default configuration (from prior use) and merge with varargin
%     varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
% end
% 
% % construct pointer to the PConn array
% CAT   = [ALLEEG.CAT];
% PConn = [CAT.PConn];
% 
% if strcmpi(typeproc,'nogui')
%     % get the config from function
%     cfg = arg_tovals(arg_report('rich',fcnHandle,[{'PConn',PConn},varargin]),false);
% else
%     % render the GUI
%     [PGh figh] = feval(['gui_' fcnName],PConn,varargin{:});
%     
%     if isempty(PGh)
%         % user chose to cancel
%         cfg = [];
%         return;
%     end
%     
%     % get the specification of the PropertyGrid
%     ps = PGh.GetPropertySpecification;
%     cfg = arg_tovals(ps,false);
% end
% 
% drawnow;
% 
% if strcmpi(typeproc,'cfg_only')
%     return;
% end
% 
% % execute the low-level function
% for cnd=1:length(ALLEEG)
%     [ALLEEG(cnd).CAT.Stats] = feval(fcnHandle,'EEG',ALLEEG(cnd),cfg);
%     
%     if ~isempty(cfg)
%         % store the configuration structure
%         ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
%     end
% end




statcondargs = '';   

% put CAT structures into array
CAT = [ALLEEG.CAT];

switch ALLEEG(1).CAT.PConn.mode 

    case 'PhaseRand'
        statcondargs = {'tail','one'};
        inputargs = {'BootstrapConnectivity',[CAT.Conn],'NullDistribution',[CAT.PConn]};
        
    otherwise
        statcondargs = {'tail','both'};
        inputargs = {'BootstrapConnectivity',[CAT.PConn]};
end
    
% render the GUI
varargin = [varargin 'statcondargs',{statcondargs}];
[PGh figh] = gui_surrogateStats([CAT.PConn],varargin{:});

if isempty(PGh)
    % user chose to cancel
    cfg = [];
    return;
end

% get the specification of the PropertyGrid
ps = PGh.GetPropertySpecification;
cfg = arg_tovals(ps,false);

drawnow

% execute the low-level function
    
% return statistics and the mean of the bootstrap estimator (if
% available)
[Stats ConnMean] = stat_surrogateStats(inputargs{:},cfg);
for cnd=1:length(ALLEEG)
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
    
    ALLEEG(cnd).CAT.Stats = Stats;
        
    if ~isempty(ConnMean(cnd))
        % replace Conn object with mean of bootstrap distribution
        % insert missing fields into new Conn object

            extrafields = setdiff(fieldnames(ALLEEG(cnd).CAT.Conn),hlp_getConnMethodNames(ALLEEG(cnd).CAT.Conn));
            for i=1:length(extrafields)
                ConnMean(cnd).(extrafields{i}) = ALLEEG(cnd).CAT.Conn.(extrafields{i});
            end 
            ALLEEG(cnd).CAT.Conn = ConnMean(cnd);
            ALLEEG(cnd).CAT.configs.TimeFreqGrid = [];
            ALLEEG(cnd).CAT.configs.BrainMovie3D = [];
    end
end



