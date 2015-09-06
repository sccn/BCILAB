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

if isfield(ALLEEG(1).CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
end

% reset the defaults if we have multiple EEG datasets
% this ensures that config selections stored in each datset (i.e. Hbase) 
% don't conflict with the allowable selection when two datasets are present
if length(ALLEEG)>1 && ~isempty(varargin)
    idx = find(ismember_bc(varargin(1:2:end),'statTest'))*2;
    if ~isempty(idx)
        for k=1:length(idx)
            varargin{idx(k)} = {}; 
        end
    end
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG},varargin]),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],ALLEEG,varargin{:});
    
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

% execute the low-level function
[Stats ConnMean] = feval(fcnHandle,'EEG',ALLEEG,cfg);

if isempty(Stats)
    % user canceled
    return;
end
    
% store statistics and (optionally) the mean of the bootstrap estimator
if length(ALLEEG)>1 && any(strcmp(cfg.statTest.arg_selection,{'Hab'}))
    % create a new dataset containing expected difference between conds
%     EEG2         = ALLEG(2);
%     EEG2.data    = -EEG2.data;  % invert data so average is mean ERP difference
%     EEG2.icaact  = -EEG2.icaact;
%     EEG_new           = pop_mergeset(ALLEEG(1),EEG2,1);
    
    
    EEG_new = ALLEEG(1);
    EEG_new.CAT = rmfield(EEG_new.CAT,'PConn');
    % compute ERP condition difference
    for fn={'data','icaact','srcpot'}
        if isfield(EEG_new,fn{1}) && ~isempty(EEG_new.(fn{1}))
            EEG_new.(fn{1}) = mean(ALLEEG(Stats.diffOrder(1)).(fn{1}),3) ...
                            - mean(ALLEEG(Stats.diffOrder(2)).(fn{1}),3);
        end
    end
    EEG_new.trials    = 1;
    EEG_new.epoch     = [];
    EEG_new.setname   = cfg.statTest.datasetOrder;
    EEG_new.condition = cfg.statTest.datasetOrder;
    EEG_new = eeg_checkset(EEG_new);
    
    EEG_new.CAT.Conn  = ConnMean;  % condition difference
    EEG_new.CAT.configs.(fcnName) = cfg;
    EEG_new.CAT.Stats = Stats;
    
    EEG_new.CAT.configs.vis_TimeFreqGrid = [];
    EEG_new.CAT.configs.vis_causalBrainMovie3D = [];
    
    % store the new EEG dataset
    [ALLEEG EEG_new] = eeg_store(ALLEEG,EEG_new,length(ALLEEG)+1);
    ALLEEG = EEG_new;
elseif length(ALLEEG)==1

    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG.CAT.configs.(fcnName) = cfg;
    end

    ALLEEG.CAT.Stats = Stats;

    if ~isempty(ConnMean)
        % replace Conn object with mean of bootstrap distribution
        % insert missing fields into new Conn object
%         extrafields = setdiff_bc(fieldnames(ALLEEG.CAT.Conn),hlp_getConnMethodNames(ALLEEG.CAT.Conn));
%         for i=1:length(extrafields)
%             ConnMean(cnd).(extrafields{i}) = ALLEEG.CAT.Conn.(extrafields{i});
%         end 
        ALLEEG.CAT.Conn = ConnMean;
        
        ALLEEG.CAT.configs.vis_TimeFreqGrid = [];
        ALLEEG.CAT.configs.vis_causalBrainMovie3D = [];
    end
end




