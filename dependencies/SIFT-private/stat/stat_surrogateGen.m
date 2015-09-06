function [PConn g] = stat_surrogateGen(varargin)
%
% Generate surrogate statistical distributions of SIFT connectivity and
% other estimators. This function can estimate bootstrap, jacknife, k-fold crossvalidation and
% phase-randomized (null) distributions.
%
%
% ===============================================
%    This function is under development and may
%    be unstable.
%    Please check sccn.uscd.edu/wiki/SIFT for 
%    updated version
% ===============================================
%
% === TODO ============================================
% - retest phaserand
% =====================================================
%
%
% Input                    Information
% ------------------------------------------------------------------------------------------------------------------------------
% ALLEEG:                  EEG data structure containing connectivity structure in EEG.CAT.PConn
%
% configs:                 Config superstructure. Contains configuration structures as fields: prepcfg,modfitcfg,conncfg. 
%                          These are returned from the respective pre_prepData(), est_fitMVAR(), and est_mvarConnectivity() 
%                          functions
%
% Optional                 Information                                                                                           
% -------------------------------------------------------------------------------------------------------------------------------
% Mode:                    Resampling modes                                                                                      
%                          Bootstrap (Efron Bootstrap resampling with replacement), Jacknife (leave-one-out cross-validation),   
%                          Crossval (k-fold cross-validation), PhaseRand (Theiler phase randomization)                                                                    
%                          Possible values: 'Bootstrap','Jacknife','SingleTrials','Crossval','PhaseRand'                      
%                          Default value  : 'Bootstrap'                                                                          
%                          Input Data Type: string                                                                               
%                                                                                                                                
%     | NumFolds:          Number of folds                                                                                       
%                          This performs k-fold resampling                                                                       
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 10                                                                                     
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
%     | NumPermutations:   Number of resamples                                                                                   
%                          This performs Theiler phase randomization                                                             
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 200                                                                                    
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
%     | NumPermutations:   Number of resamples                                                                                   
%                          This performs efron bootstrapping                                                                     
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 200                                                                                    
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
% AutoSave:                Autosave distributions                                                                                
%                          This will periodically save the computation to a mat file. In case of system crash, bootstrapping     
%                          can be resumed from this file                                                                         
%                          Input Range  : Unrestricted                                                                           
%                          Default value: 0                                                                                      
%                          Input Data Type: boolean                                                                              
%                                                                                                                                
%     | FileNamePrefix:    Prefix (including optional path) for autosave file                                                    
%                          Data is saved in <FileNamePrefix>.part                                                                
%                          Possible values: Unrestricted                                                                         
%                          Default value  : 'SIFT_boostrap'                                                                      
%                          Input Data Type: string                                                                               
%                                                                                                                                
%     | AutoSaveFrequency: Fractional increment between saves                                                                    
%                          For example, 0.25 = quarterly saves                                                                   
%                          Input Range  : [2.2204e-16           1]                                                               
%                          Default value: 0.25                                                                                   
%                          Input Data Type: real number (double)                                                                 
%                                                                                                                                
% ConnectivityMethods:     Connectivity estimator(s) to bootstrap                                                                
%                          All connectivity estimators must be named exactly as they appear in EEG.CAT.Conn                      
%                          Possible values: ''                                                                                   
%                          Default value  : 'n/a'                                                                                
%                          Input Data Type: boolean                                                                              
%                                                                                                                                
% VerbosityLevel:          Verbosity level. 0 = no output, 1 = text, 2 = graphical                                               
%                          Possible values: 0,1,2                                                                                
%                          Default value  : 2                                                                                    
%                          Input Data Type: real number (double) 
%
% Output                 Information                                                                                           
% -------------------------------------------------------------------------------------------------------------------------------
% PConn:                 Surrogate connectivity structure. Same format as EEG.CAT.Conn but with resamples stored in the last 
%                        dimension of each connectivity matrix (e.g. [M x M x Time x Freq x Resamples])   
%
% See Also: stat_analyticStats(), statcond(), stat_surrogateStats()
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 5.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD.
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



% extract some stuff from inputs for arg defaults
ALLEEG = arg_extract(varargin,{'EEG','ALLEEG'},1);

if length(ALLEEG)>1
    error('Surrogate distributions must be computed for each dataset individually');
end

% defaults
mvarConnectivityDef = {};
modelingDef         = hlp_getModelingApproaches('defNameOnly');

if ~isempty(ALLEEG)
    
    % if model was previously fit, set default modeling
    % approach from the model
    if isempty(hlp_checkeegset(ALLEEG,{'model'}))
        modelingDef = ALLEEG.CAT.MODEL.modelapproach;
    end
    
    % get the name of of the modeling approach function
    modFcnName  = hlp_getModelingApproaches('mfileNameOnly',modelingDef);

    % get the last-used configuration structures...
    configs = ALLEEG.CAT.configs;
    
    % ... and update defaults (if possible)
    if ~isempty(configs.est_mvarConnectivity)
        mvarConnectivityDef = {configs.est_mvarConnectivity};
    end
    
    if ~isempty(configs.(modFcnName))
        modelingDef = {modelingDef configs.(modFcnName)};
    end
    
end

clear ALLEEG Conn;

g = arg_define([0 Inf],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEG structure.'), ...
    arg_subswitch({'modelingApproach','ModelingApproach'},modelingDef, ...
        hlp_getModelingApproaches, ...
        'Select a modeling approach and define parameters.', ...
        'cat','Modeling Pipeline', ...
        'suppress',{'verb','timer'}), ...
    arg_sub({'connectivityModeling','ConnectivityEstimation'},mvarConnectivityDef,@est_mvarConnectivity,'Connectivity estimation options','suppress',{'verb'},'cat','Modeling Pipeline'), ...
    arg_subswitch({'mode','Mode'},'Bootstrap', ...
        { ...
        'Bootstrap' { ...
        arg({'nperms','NumPermutations'},200,[],'Number of resamples. This performs efron bootstrapping'), ...
        arg({'saveTrialIdx','SaveTrialIndices'},true,[],'Save trial indices for each resample iteration.'), ...
        }, ...
        'Jacknife', { ...
        arg({'saveTrialIdx','SaveTrialIndices'},true,[],'Save trial indices for each resample iteration.'), ...
        }, ...
        'SingleTrials', { ...
        arg({'saveTrialIdx','SaveTrialIndices'},true,[],'Save trial indices for each resample iteration.'), ...
        }, ...
        'Crossval' { ...
        arg({'nfolds','NumFolds'},10,[],'Number of folds. This performs k-fold resampling'), ...
        arg({'saveTrialIdx','SaveTrialIndices'},true,[],'Save trial indices for each resample iteration.'), ...
        }, ...
        'PhaseRand', { ...
        arg({'nperms','NumPermutations'},200,[],'Number of resamples. This performs Theiler phase randomization'), ...
        }, ...
        'CustomSchedule', { ...
        arg({'resampleTrialIdx','ResamplingSchedule'},[],[],'Resampling schedule matrix [num_iterations x sample_size]. Each row is a vector of trial indices which constitute the sample set for a single resampling iteration.','shape','matrix'), ...
        }}, ...
        'Resampling modes. Bootstrap (Efron Bootstrap resampling with replacement), Jacknife (leave-one-out cross-validation), Crossval (k-fold cross-validation), PhaseRand (Theiler phase randomization)', ...
        'cat','Resampling Options' ...
        ), ...
    arg_subtoggle({'autosave','AutoSave'},[], ...
        { ...
        arg({'savefname','FileNamePrefix'},'SIFT_boostrap','','Prefix (including optional path) for autosave file. Data is saved in <FileNamePrefix>.part'), ...
        arg({'savefrq','AutoSaveFrequency'},0.25,[eps 1-eps],'Fractional increment between saves. For example, 0.25 = quarterly saves'), ...
        }, 'Autosave distributions. This will periodically save the computation to a mat file. In case of system crash, bootstrapping can be resumed from this file', ...
        'cat','Resampling Options'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );


arg_toworkspace(g);
ALLEEG = EEG; clear EEG;

% check the dataset
chk = hlp_checkeegset(ALLEEG,{'cat'});
if ~isempty(chk)
    error(chk{1}); 
end

if (ALLEEG.CAT.trials==1 || isempty(ALLEEG.CAT.trials)) ...
        && any(ismember_bc(lower(mode.arg_selection),{'bootstrap','jacknife','inversejacknife','crossval'}))
    error(['Unable to compute bootstrap distributions for a single trial. ' char(10) ...
           'You must use either Analytic Statistics or Phase Randomization surrogate option']);
end
    
% Perform boostrapping
switch mode.arg_selection
    case 'Bootstrap'
        nperms = mode.nperms;
        cv.training = @(x)(randsample(ALLEEG.CAT.trials, ALLEEG.CAT.trials, true));
    case 'Jacknife'
        nperms = ALLEEG.CAT.trials;
        cv = cvpartition(ALLEEG.CAT.trials,'leaveout');
    case 'SingleTrials'
        nperms = ALLEEG.CAT.trials;
%         cv = cvpartition(ALLEEG.CAT.trials,'leaveout');
    case 'Crossval'
        nperms = mode.nfolds;
        cv = cvpartition(ALLEEG.CAT.trials,'kfold',nperms);
    case 'PhaseRand'
        nperms = mode.nperms;
        mode.saveTrialIdx = false;
    case 'CustomSchedule'
        nperms = size(mode.resampleTrialIdx,1);
        mode.saveTrialIdx = true;
end

% If necessary, init matrix for storing resampling trial indices
if mode.saveTrialIdx
    switch mode.arg_selection
        case 'SingleTrials'
            resampleTrialIdx = zeros(nperms,nperms,'int32');
        case {'Bootstrap','Jacknife','Crossval'}
            resampleTrialIdx = zeros(nperms,length(cv.training(1)),'int32');
        case 'CustomSchedule'
            resampleTrialIdx = mode.resampleTrialIdx;
        otherwise
            resampleTrialIdx = [];
    end
else
    resampleTrialIdx = [];
end

% Initiate progress bar
if verb==1
    progress('init',sprintf('Resampling (%d/%d)...',0,nperms));
    tic
    t0 = toc;
elseif verb==2
    waitbarTitle = sprintf('%s Resampling (x%d)...',mode.arg_selection, nperms);
    multiWaitbar(waitbarTitle, ...
                 'Color', [1.0 0.4 0.0], ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end

if autosave.arg_selection
    savefrq = round(autosave.savefrq*nperms);
end


% get the m-file name of the function implementing the modeling approach
modelingFuncName = hlp_getModelingApproaches('mfileNameOnly',g.modelingApproach.arg_selection);

% main bootstrap loop
% -------------------------------------------------------------------------
for perm=1:nperms
    
    % restore the original dataset
    EEG = ALLEEG;
            
    % select the trials for this resample
    switch mode.arg_selection
        case 'SingleTrials'
            % fit model to each single trial
            rsmpIdx = perm; %cv.test(perm);
        case 'PhaseRand'
            % phase randomization
            rsmpIdx = [];
        case 'CustomSchedule'
            % use predetermined trial indices
            rsmpIdx = mode.resampleTrialIdx(perm,:);
        otherwise
            % all other bootstrapping/resampling modes
            rsmpIdx = cv.training(perm);
    end
    
    if strcmp(mode.arg_selection,'PhaseRand')
        
        % Create surrogate data by randomizing phases.
        % Multiply each fourier amplitude by e^{iw)
        % where w is a random phase chosen in [0 2pi]
        % (c.f. Theiler, et al 1997)
        % To generate hermitian phase distributions, we extract
        % the phase of a random matrix. This ensures the surrogate
        % spectrum is conjugate symmetric
        EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1 3]);
        [npnts nchs ntr] = size(EEG.CAT.srcdata);
        for tr=1:ntr
            EEG.CAT.srcdata(:,:,tr) = ...
                ifft(abs(fft(EEG.CAT.srcdata(:,:,tr))) ...
                .* exp(1i*angle(fft(rand(npnts,nchs)))), ...
                'symmetric');
        end
        EEG.CAT.srcdata = ipermute(EEG.CAT.srcdata,[2 1 3]);
        
    else
        % select subset of trials
        EEG.CAT.srcdata = EEG.CAT.srcdata(:,:,rsmpIdx);
        EEG.CAT.trials  = length(rsmpIdx);
    end
    
    if mode.saveTrialIdx
        resampleTrialIdx(perm,:) = rsmpIdx;
    end
        
    % Fit model and estimate connectivity for this permutation
    % ---------------------------------------------------------------------
   
    % fit model
    MODEL = feval(modelingFuncName,'EEG',EEG,g.modelingApproach,'verb',0);
    
    % calculate connectivity
    C = est_mvarConnectivity('EEG',EEG,'MODEL',MODEL,g.connectivityModeling,'verb',0);

    % append causal estimates to the last dimension of respective arrays in Conn
    connfields = hlp_getConnMethodNames(C);
    if perm==1
        for m=1:length(connfields)
            if ~isempty(C.(connfields{m}))
                % initialize matrices
                % permutations will be temporarily stored in first
                % dimension
                PConn.(connfields{m}) = zeros([nperms size(C.(connfields{m}))],'single');
            end
        end
    end
    for m=1:length(connfields)
        if ~isempty(C.(connfields{m}))
            % store permutation
            PConn.(connfields{m})(perm,:,:,:,:) = single(C.(connfields{m}));
        end
    end

    % autosave checkpoint
    if autosave.arg_selection && ~isempty(autosave.savefname) && ~mod(perm,savefrq)
        save(sprintf('%s.part',autosave.savefname),'PConn','-v7.3');
    end
    
    % ---------------------------------------------------------------------
    
    % update progress bar
    if verb==1
        % text progress bar
        te = toc-t0;
        tt = ceil(te*(nperms-perm)/perm);
        progress(perm/nperms, ...
            sprintf('Resampling [%d/%d] (EL: %0.1fm / ETA: %0.1fm)', ...
            perm,nperms,ceil(te)/60,tt/60));
    elseif verb==2
        % graphical waitbar
        cancel = multiWaitbar(waitbarTitle,perm/nperms);
        if cancel
            if strcmpi('yes',questdlg2( ...
                            'Are you sure you want to cancel?', ...
                            'Resampling','Yes','No','No'));
                PConn = [];
                multiWaitbar(waitbarTitle,'Close');
                return;
            else
                multiWaitbar(waitbarTitle,'ResetCancel',true);
            end
        end 
    end
    
end

% construct PConn object
PConn.winCenterTimes    = C.winCenterTimes;
PConn.erWinCenterTimes  = C.erWinCenterTimes;
PConn.freqs             = C.freqs;
PConn.mode              = mode.arg_selection;
PConn.resampleTrialIdx  = resampleTrialIdx;
for m=1:length(connfields)
    sz = size(PConn.(connfields{m}));
    % put permutations in last dimension
    PConn.(connfields{m}) = permute(PConn.(connfields{m}),circshift((1:length(sz))',-1));
    
%     % check dimensions
%     szp = size(PConn.(connfields{m}));
%     szs = size(EEG.CAT.Conn.(connfields{m}));
%     [dummy dimidx] = setdiff_bc(szp(1:end-1),szs);
%     if ~isempty(dimidx)
%         % a singleton dimension was squeezed out, restore it
%         PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),dimidx+1);
%     end
    
    % insert singleton dimensions if necessary
    if length(PConn.winCenterTimes)==1
        PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),4);
    end
    if length(PConn.freqs)==1
        PConn.(connfields{m}) = hlp_insertSingletonDim(PConn.(connfields{m}),3);
    end
end

    
    
    

% clean up
if verb==1
    fprintf('\n')
elseif verb==2
    multiWaitbar(waitbarTitle,'Close');
end

if autosave.arg_selection ...
   && ~isempty(autosave.savefname) ...
   && exist(sprintf('%s.part',autosave.savefname),'file')
   % delete autosave file
    delete(sprintf('%s.part',autosave.savefname));
end
