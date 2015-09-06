

function res = vis_eigmodes(varargin)
%
% Visualize the results of an eigenmode decomposition of a VAR[p] model as 
% computed using est_eigmode().
% This plots the characteristic frequencies and damping times of the
% eigenmodes, optionally sorted by dynamical importance [1].
%
% TODO: this function is incomplete. Plotting functionality not yet
% implemented in this function.
%
% Inputs:
%
%   EEG:        EEGLAB dataset containing CAT.EIGMODE substructure
%   
% Optional:
% 
%     TopPercentEigenmodes: Top percent eigenmodes (0-100)                                                                        
%                           The percent of eigmodes with highest dynamical importance to plot                                     
%                           Input Range  : [0  100]                                                                               
%                           Default value: 100                                                                                    
%                           Input Data Type: real number (double)                                                                 
% 
%     ConvertPeriodToHz:    Convert res.period to Hertz                                                                           
%                           Input Data Type: boolean                                                                              
% 
%     sortEigmodes:         Sort eigmodes by dynamical importance                                                                 
%                           If true, sort for each window independently.If false, and topPercent < 100, sort each window,         
%                           determine the topPercent most frequently occuring eigenmodes (mode across windows),and return these   
%                           eigenmodes, unsorted across windows                                                                   
%                           Input Data Type: boolean                                                                              
% 
%     plotResults:          Plot results                                                                                          
%                           Input Data Type: boolean 
% 
% Output:
% 
%   res:        results structure
%
% See Also: est_eigenmode()
%
% References:
% 
% [1] Neumaier A, Schneider T (2001) Estimation of parameters and 
%   eigenmodes of multivariate autoregressive models. ACM Transactions on 
%   Mathematical Software (TOMS) 27:57
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



    g = arg_define([0 1],varargin, ...
            arg_norep('EEG',mandatory), ...
            arg({'topPercent','TopPercentEigenmodes'},100,[0 100],'Top percent eigenmodes (0-100). The percent of eigmodes with highest dynamical importance to plot'), ...
            arg({'convertPeriodToHz','ConvertPeriodToHz'},true,[],'Convert res.period to Hertz.'), ...
            arg({'sortEigmodes'},true, [], ['Sort eigmodes by dynamical importance. ' ...
                                            'If true, sort for each window independently.' ...
                                            'If false, and topPercent < 100, sort each window, determine the topPercent most frequently occuring eigenmodes (mode across windows),' ...
                                            'and return these eigenmodes, unsorted across windows']), ...
            arg({'plotResults'},false,[],'Plot results') ...
            );
    
    % commit some variables to the workspace
    [data g] = hlp_splitstruct(g,{'EEG'});        
    arg_toworkspace(data);
    clear data;
    
    
    [nch nmodes] = size(EEG.CAT.EIGMODE.modes{1});
    nwins        = length(EEG.CAT.EIGMODE.modes);
    
    % convert cell data to matrix format
    res.eigmodes    = cell2mat(EEG.CAT.EIGMODE.modes);
    res.eigmodes    = reshape(res.eigmodes,[nch nmodes nwins]);
    
    if ~isempty(EEG.CAT.EIGMODE.modeconf{1})
        res.eigmodesConf = cell2mat(EEG.CAT.EIGMODE.modeconf);
        res.eigmodesConf = reshape(res.eigmodesConf,[nmodes nwins]);
    else
        res.eigmodesConf = zeros(nmodes, nwins);
    end
        
    res.period = cell2mat(EEG.CAT.EIGMODE.period);
    res.period = reshape(res.period,[2 nmodes nwins]);
    res.periodConf = squeeze(res.period(2,:,:));
    res.period = squeeze(res.period(1,:,:));
    
    res.dampingTime = cell2mat(EEG.CAT.EIGMODE.dampingTime);
    res.dampingTime = reshape(res.dampingTime,[2 nmodes nwins]);
    res.dampingTimeConf = squeeze(res.dampingTime(2,:,:));
    res.dampingTime = squeeze(res.dampingTime(1,:,:));
    
    res.eigenvalues      = cell2mat(EEG.CAT.EIGMODE.lambda);
    
    % sort the eigmodes by dynamical importance
    res.exctn = cell2mat(EEG.CAT.EIGMODE.exctn);
    res.exctn = reshape(res.exctn,[nmodes nwins]);

    % sort eigmodes by dynamical importance
    if g.sortEigmodes
        [dummy res.topEigmodesIndex] = sort(res.exctn,1,'descend');
    elseif g.topPercent < 100
        % determine which modes occur most frequently across windows
        % and sort rows -- for all windows jointly -- in this order
        [dummy topmodes] = sort(res.exctn,1,'descend');
        topmodes = mode(topmodes,2);
        
        % remove duplicate topmodes
        for i=1:length(topmodes)
            val = topmodes(i);
            if ~isnan(val)
                topmodes(topmodes==topmodes(i))=NaN; % mark all duplicates
                topmodes(i)=val;
            end
        end
        topmodes(isnan(topmodes))=[];
        
%         [b m n] = unique_bc(topmodes);
%         mostFreqModes = topmodes(sort(m));
        res.topEigmodesIndex = repmat(topmodes,1,nwins);
    else
        res.topEigmodesIndex = repmat((1:nmodes)',1,nwins);
    end
    
    % extract the top percent of most important eigenmodes
    numTopEigmodes = round((g.topPercent/100)*size(res.topEigmodesIndex,1));
    res.topEigmodesIndex = res.topEigmodesIndex(1:numTopEigmodes,:);
    
    % extract the data corresponding to the top eigenmodes
    colshift = nmodes*(0:nwins-1);
    eigmodeLinearIndex = res.topEigmodesIndex + repmat(colshift,numTopEigmodes,1);
    res.exctn       = res.exctn(eigmodeLinearIndex);
    res.eigmodesConf = res.eigmodesConf(eigmodeLinearIndex);
    res.period       = res.period(eigmodeLinearIndex);
    res.periodConf   = res.periodConf(eigmodeLinearIndex);
    res.dampingTime  = res.dampingTime(eigmodeLinearIndex);
    res.dampingTimeConf = res.dampingTimeConf(eigmodeLinearIndex);
    res.eigenvalues     = res.eigenvalues(eigmodeLinearIndex);
    tmp = zeros(nch,numTopEigmodes,nwins);
    for ch=1:nch
        tmpeigmode  = squeeze(res.eigmodes(ch,:,:));
        tmp(ch,:,:) = tmpeigmode(eigmodeLinearIndex);
    end
    res.eigmodes = tmp;
    
    res.convertPeriodToHz = g.convertPeriodToHz;
    if g.convertPeriodToHz
        res.period = 1./res.period;
        res.periodConf = 1./res.periodConf;
    end
    
