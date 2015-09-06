function [FIT MCMC_LastState]=stat_ubsplsm(varargin)
% univariate b-spline smoothing with fPCA

% Y = cell array of [numChans x numChans x T] connectivities for each subject
% K = number of knots
% Q = number of fpca basis functions

% output:
% fit_distrib: cell array of [numChans x numChans x T x iterations]
%              posterior distribution for each subject

% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
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

arg_define([0 1],varargin, ...
    arg_norep({'Y','TSData'},mandatory,[],sprintf(['Cell array of data to smooth.\n' ...
              'Generally, Y is a cell array where Y{i} is the T x 1 vector of time-varying (or freq-varying) connectivity for the kth channel pair of the sth subject.\n' ...
              'e.g. Y = {s1(1,1) s1(1,2) s1(1,3) ... s1(N,1) s1(N,2) s1(N,3) ... \n' ...
              '          s2(1,1) s2(1,2) s2(1,3) ... } \n'])), ...
    arg({'smoothingLayout','MatrixElementsToSmooth'},{'diagonals','off-diagonals'},{'diagonals','off-diagonals'},'Which parts of the matrix to smooth. Diagonals (e.g. auto-connectivity) and off-diagonals (e.g. cross-connectivity) will be smoothed separately','type','logical'), ...
    arg_sub({'mcmc_opts','MCMC'},{}, ...
    {...
        arg_sub({'mcmc_init','InitMCMC'},{},@stat_ubsplsm_init,'Initialize MCMC','suppress','verb') ...
        arg_sub({'mcmc_run','RunMCMC'},{},@stat_ubsplsm_mcmc,'MCMC runtime options','suppress','verb') ...
    },'Perform Markov Chain Monte Carlo estimation'), ...
    arg_norep({'MCMC_InitState'},struct([]),[],'Object containing initial state of Gibbs sampler. If supplied, this overrides mcmc_init'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );
    
% initialize some vars
init_mcmc = isempty(g.MCMC_InitState);    
nsubj=length(Y);

% determine the number of variables for each subject
m_i=zeros(nsubj,1);
for i=1:nsubj
    m_i(i)=size(Y{i},1);
end

% determining smoothing layout
smoothDiags     = ismember_bc('diagonals',g.smoothingLayout);
smoothOffDiags  = ismember_bc('off-diagonals',g.smoothingLayout);
    
%% get smoothed time-varying connectivity coefficients

% Diagonals and off-diagonals often have different variances, so we smooth
% them separately

% smooth diagonals
% -------------------------------------------------------------------------
if smoothDiags
    if verb==2
        multiWaitbar('Smoothing Diagonals','Reset', ...
                     hlp_getNextUniqueColor('reset'));
    end

    % extract self-connectivity for each subject
    CPairs=cell(sum(m_i),1);
    ind=0;
    for i=1:nsubj
        for j=1:m_i(i)
            ind=ind+1;
            CPairs{ind}=squeeze(Y{i}(j,j,:));
        end
    end

    % initialize MCMC
    if init_mcmc
        MCMC_InitState = stat_ubsplsm_init('Y',mcmc_opts.mcmc_init,'verb',verb);
    end
    
    % perform the smoothing
    [fit_diag MCMC_LastState.diag]=stat_ubsplsm_mcmc('Y',CPairs,mcmc_opts.mcmc_run, ...
                                        'verb',verb,'MCMC_InitState',MCMC_InitState.diag);
end


% smooth off-diagonals
% -------------------------------------------------------------------------
if smoothOffDiags
    if verb
        multiWaitbar('Smoothing Off-Diagonals',1/3);
    end

    % count the total number of pairs
    m_offdiag=sum(m_i.^2-m_i);

    % extract cross-connectivity off-diagonal
    % elements for each subject
    CPairs=cell(m_offdiag,1);
    ind=0;
    for i=1:nsubj
        for j1=1:m_i(i)
            for j2=1:m_i(i)
                if j1~=j2
                    ind=ind+1;
                    CPairs{ind}=squeeze(Y{i}(j1,j2,:));
                end
            end
        end
    end

    % initialize MCMC
    if init_mcmc
        MCMC_InitState = stat_ubsplsm_init('Y',mcmc_opts.mcmc_init,'verb',verb);
    end
    
    % perform the smoothing
    [fit_offdiag MCMC_LastState.offdiag]=stat_ubsplsm_mcmc('Y',CPairs,mcmc_opts.mcmc_run, ...
                                         'verb',verb,'MCMC_InitState',MCMC_InitState.offdiag);
end


% store data in [nchs x nchs x time x distribution]
if verb
    multiwaitbar('Creating final data matrices...',2/3);
end
FIT = cell(1,nsubj);
ind=0;
ind_diag = 0;
for i=1:nsubj
    FIT{i} = zeros(m_i(i),m_i(i),size(fit_offdiag,1),niterToKeep);
    
    for j1=1:m_i(i)
        for j2=1:m_i(i)
            if j1~=j2
                ind=ind+1;
                FIT{i}(j1,j2,:,:)=squeeze(fit_offdiag(:,ind,:));
            elseif niters~=0
                ind_diag = ind_diag+1;
                FIT{i}(j1,j2,:,:)=squeeze(fit_diag(:,ind_diag,:));
            end
        end
    end
end

if verb
    multiWaitbar('CloseAll');
    pause(0.1);
end

