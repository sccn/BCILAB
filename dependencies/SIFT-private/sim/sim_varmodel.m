function [data truth A C] = sim_varmodel(varargin)
% simulate a vector autoregressive model
%
% This returns an EEG dataset or [nchs x npnts x ntr] matrix of simulated 
% data
%
% The function also can return the 'ground truth' in an EEGLAB dataset
%
% See Also: sim_simulateSources()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2011-2013, SCCN/INC, UCSD.
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

g = arg_define([0 Inf],varargin,...
    arg_subswitch({'sim','Simulation'},hlp_getSimExamples('defaultNameOnly'), ...
    hlp_getSimExamples, ...
    {'Select a simulation.',hlp_microcache('sift_domain',@hlp_buildSimHelpText)}), ...
    arg_sub({'simParams','SimParams'},{},...
    { ...
    arg({'srate','SamplingRate'},100,[1 Inf],'Process sampling rate (Hz)'), ...
    arg({'Nl','TrialLength'},5,[1 Inf],'Trial length (sec)'), ...
    arg({'Nr','NumTrials'},100,[1 Inf],'Number of trials (realizations)'), ...
    arg({'ndisc','BurnInSamples'},1000,[0 Inf],'Burn-in samples. Number of initial simulated samples to discard') ...
    arg({'checkStability','CheckStability'},true,[],'Check whether process is stable'), ...
    },'Simulation Parameters'), ...
    arg_sub({'genParams','DataGenParams'},{'A',[]},@sim_genTVARdata,'Data generation options','suppress',{'ndisc','Nl','Nr'}), ...
    arg_subtoggle({'makeEEGset','BuildEEGLABStructure'},'on',...
    { ...
        arg_subtoggle({'exportGroundTruth','ExportGroundTruth'},'off',...
            { ...
                arg({'winlen','WindowLength'},0.5,[eps Inf],'Effective Sliding window length (sec)'), ...
                arg({'winstep','WindowStepSize'},0.03,[eps Inf],'Effective Window step size (sec)'), ...
            },'Export ground truth model. This will return an EEGLAB dataset containing the Model in SIFT format'), ...
        arg({'sourceNames','SourceNames'},[],[],'Assign labels to sources. Default is numeric index.','type','cellstr','shape','row'), ...
        arg({'setname','SetName'},'','','Dataset name. If empty, simulation name is used','type','char') ...
    },'Build EEGLAB datastructure. This will return an EEGLAB datastructure containing the simulated data. Otherwise, the raw data is returned','cat','OutputFormat'), ...
    arg({'plotData','PlotData'},true,[],'Plot simulated data','cat','Visualization'), ...
    arg({'plotGraphicalModel','PlotGraphicalModel'},true,[],'Plot graphical model (if available). The graphical model image should be rendered in a jpg and stored in the <sift-root>/resources/sim/ folder and the name should be <sim_name>.jpg where <sim_name> is the name of the simulation m-file following the prefix ''sim_ex_''','cat','Visualization'), ...
    arg({'verb','VerbosityLevel'},int32(2),{int32(0) int32(1) int32(2)},'Verbose','type','int32'));

% do some input checking
if isfield(g.sim,'morder')
    ModelOrder = g.sim.morder;
elseif isfield(g.sim,'rndopts')
    ModelOrder = g.sim.rndopts.morder;
end
if isempty(ModelOrder)
    error('SIFT:sim_varmodel:badParam','ModelOrder must be manually specified');
end
% convert from seconds to samples
g.simParams.Nl = round(g.simParams.Nl*g.simParams.srate);
% initialize outputs
[data truth]  = deal([]);
% Intialize progress bar
curstep  = 0;
numsteps = 3+g.makeEEGset.arg_selection;
waitbarTitle = sprintf('Simulating %s model...',g.sim.arg_selection);
createMasterWaitbar();

% create prototype VAR structure
% -------------------------------------------------------------------------
if isfield(g.sim,'expr')
    wbt = 'Translating system equations....';
    createWaitbar(wbt);
    Aproto = sim_genVARModelFromEq(g.sim.expr,ModelOrder);
    updateWaitbar(wbt,1);
    if ~updateMasterWaitbar(),  return; end
elseif isfield(g.sim,'rndopts')
    wbt = 'Generating random model coefficients....';
    createWaitbar(wbt);
    Aproto = sim_genRndVARcoeffs(g.sim.rndopts);
    updateWaitbar(wbt,1);
    if ~updateMasterWaitbar(),  return; end
elseif isempty(g.sim.expr)
    error('sim_varmodel:badexpr','System equations expression was empty!');
end
% generate the VAR coefficients
% -------------------------------------------------------------------------
wbt = 'Building model....';
createWaitbar(wbt);
[A stable] = sim_genTVARcoeffs('Aproto',Aproto, 'ModelOrder',ModelOrder,     ...
                               'Nl',g.simParams.Nl, 'ndisc',g.simParams.ndisc, ...
                               'checkStability',g.simParams.checkStability,    ...
                               'verb',g.verb);
% stability warning
if ~g.verb && ~all(stable)
    if length(stable)>1
        fprintf('WARNING: System is unstable!\n');
    else
        fprintf('WARNING: System is unstable at sample(s) %s!\n',hlp_tostring(find(~stable)));
    end
end
updateWaitbar(wbt,1);
if ~updateMasterWaitbar(), return; end

% generate data from the VAR model
% -------------------------------------------------------------------------
wbt = 'Simulating data from model....';
createWaitbar(wbt);
[data C mu] = sim_genTVARdata('A',A, g.genParams, 'Nl',g.simParams.Nl,   ...
                           'Nr',g.simParams.Nr, 'ndisc',g.simParams.ndisc);

updateWaitbar(wbt,1);
if ~updateMasterWaitbar(), return; end

% construct EEG dataset
% -------------------------------------------------------------------------
if g.makeEEGset.arg_selection
    wbt = 'Constructing EEG dataset....';
    createWaitbar(wbt);
    
    if isempty(g.makeEEGset.setname)
        if isfield(g.sim,'savesim')
            g.makeEEGset.setname = g.sim.savesim.simName;
        else
            g.makeEEGset.setname = g.sim.arg_selection;
        end
    end
    
    M = size(data,1);
    
    EEG             = eeg_emptyset;
    EEG.data        = data;
    EEG.icaweights  = [];
    EEG.icasphere   = [];
    EEG.icawinv     = [];
    EEG.icaact  	= [];
    EEG.srate       = g.simParams.srate;
    EEG.times       = ((0:(g.simParams.Nl-1))/g.simParams.srate)*1000; % ms
    EEG.pnts        = g.simParams.Nl;
    EEG.trials      = g.simParams.Nr;
    EEG.xmin        = EEG.times(1);
    EEG.xmax        = EEG.times(end)/1000;  % sec
    EEG.nbchan      = M;
    EEG.setname     = g.makeEEGset.setname;
    EEG.condition   = EEG.setname;
    % validate the eeg dataset
    EEG = eeg_checkset(EEG);
    % create SIFT datastructure
    if isempty(g.makeEEGset.sourceNames)
        g.makeEEGset.sourceNames = cellstr(num2str((1:EEG.nbchan)'))';
    end
    EEG.CAT         = hlp_sift_emptyset(...
                    'srcdata',  EEG.data,     ...
                    'nbchan',   EEG.nbchan,   ...
                    'pnts',     EEG.pnts,     ...
                    'trials',   EEG.trials,   ...
                    'times',    EEG.times,    ...
                    'curComps', 1:EEG.nbchan, ...
                    'curComponentNames',g.makeEEGset.sourceNames);
    % return EEG data
    data = EEG;
    
    updateWaitbar(wbt,1);
    
    
    if g.makeEEGset.exportGroundTruth.arg_selection
        % construct ground truth dataset as well
        wbt = 'Constructing Ground Truth EEG dataset....';
        createWaitbar(wbt);
        truth = makeOracle(EEG,A,C,mu,ModelOrder, ...
                            g.makeEEGset.exportGroundTruth.winstep, ...
                            g.makeEEGset.exportGroundTruth.winlen);
        updateWaitbar(wbt,1);
    end
    
    if ~updateMasterWaitbar(), return; end
end

% plot the data
% -------------------------------------------------------------------------
if g.plotData
    if isstruct(data)
        pop_eegplot(data);
    else
        eegplot(data,'srate',g.simParams.srate);
    end
end

% plot the graphical model for this example
% -------------------------------------------------------------------------
if g.plotGraphicalModel
    % obtain the resource name from the suffix of the sim m-file
    simExs = hlp_getSimExamples;
    simIdx = cellfun(@(x) strcmp(x{1},g.sim.arg_selection),simExs);
    simFcnName = func2str(simExs{simIdx}{2});
    simRscName = strrep(simFcnName,'sim_ex_','');
    try
        hlp_viewGraphicsResource(['sim/' simRscName '.jpg']);
    catch err
        fprintf('Failed to plot graphical model. Got error instead: %s\n',err.message);
    end
end

% save the simulation as an m-file
% -------------------------------------------------------------------------
if isfield(g.sim,'savesim') && g.sim.savesim.arg_selection
    writeSim(g.sim.expr,g.sim.morder, ...
             g.sim.savesim.simName,   ...
             g.sim.savesim.simDescrip,...
             g.sim.savesim.simAuthor, ...
             g.sim.savesim.simRef,    ...
             g.verb);
end

% finally, cleanup
deleteAllWaitbars();


% helper functions
% -------------------------------------------------------------------------

    function createMasterWaitbar()
        if g.verb==1
            fprintf('---------------\n');
            fprintf('%s\n',waitbarTitle);
            fprintf('---------------\n');
        elseif g.verb==2
            % create waitbar
            multiWaitbar(waitbarTitle,'Reset');
            multiWaitbar(waitbarTitle,'ResetCancel',true);
            multiWaitbar(waitbarTitle, ...
                'Color', hlp_getNextUniqueColor('reset'), ...
                'CanCancel','on', ...
                'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
        end
    end

    function nocancel=updateMasterWaitbar()
        nocancel = true;
        if g.verb==1
            fprintf('---------------\n');
        elseif g.verb==2
            curstep = curstep + 1;
            drawnow;
            cancel = multiWaitbar(waitbarTitle,curstep/numsteps);
            if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
                deleteAllWaitbars()
                nocancel = false;
                return;
            end
        end
    end

    % createWaitbar
    function createWaitbar(titleString)
        if g.verb==1
            fprintf('%s\n',titleString);
        elseif g.verb==2
            multiWaitbar(titleString, ...
                'Color', hlp_getNextUniqueColor(), ...
                'CanCancel','off');
        end
    end

    % deleteWaitbar
    function deleteWaitbar(titleString)
        if g.verb==1
            fprintf('\n');
        elseif g.verb==2
            multiWaitbar(titleString,'Close');
        end
    end

    % updateWaitbar
    function updateWaitbar(titleString,fracDone)
        if g.verb==1
            fprintf('%0.2f%% ',fracDone*100);
            if (fracDone == 1)
                fprintf('\n'); end
        elseif g.verb==2
            multiWaitbar(titleString,fracDone);
        end
    end

    % deleteWaitbar
    function deleteAllWaitbars()
        if g.verb==2
            multiWaitbar('CloseAll');
        end
    end

    % makeOracle. This function returns a "oracle" EEG dataset
    % containing a ground truth model
    % A is the matrix or cell array containing VAR coefficients
    % C is the noise covariance matrix
    function EEGtrue = makeOracle(EEGsim,A,C,mu,ModelOrder,winStep,winLen)
        
        % first do some bookeeping
        if rem(winStep,1/EEGsim.srate)
            if winStep<1/EEGsim.srate
                winStep = 1/EEGsim.srate;
            else
                % adjust step size to nearest multiple of sampling interval
                winStep = winStep-rem(winStep,1/EEGsim.srate);
            end
            if g.verb,
                fprintf('Adjusting window step size to nearest multiple of sampling interval\n');
                fprintf('\tstep size is now %0.5g sec\n',winStep);
            end
        end
        if rem(winLen,1/EEGsim.srate)
            if winLen<1/EEGsim.srate
                winLen = 1/EEGsim.srate;
            else
                % adjust window length to nearest multiple of sampling interval
                winLen = winLen-rem(winLen,1/EEGsim.srate);
            end
            if g.verb,
                fprintf('Adjusting window length to nearest multiple of sampling interval\n');
                fprintf('\twindow length is now %0.5g sec\n',winLen);
            end
        end
        if winLen > EEGsim.xmax-EEGsim.xmin
            winLen = EEGsim.xmax-EEGsim.xmin;
            if g.verb,
                fprintf('Window length exceeds trial length. Adjusting window length to match trial length\n');
                fprintf('\twindow length is now %0.5g sec\n',winLen);
            end
        end
        tidx = [1 EEGsim.pnts];

        winLenPnts  = round(winLen*EEGsim.srate); % window size in points
        winStepPnts = round(winStep*EEGsim.srate);% step size in points

        % starting point of each window (points)
        winStartIdx  = tidx(1):winStepPnts:(tidx(2)-winLenPnts)+1;
        winStartTimes = (winStartIdx-1)/EEGsim.srate;
        
        % construct true AR model object
        
%         wintimes = winStartTimes+(winLen)/2;
%         idx = getindex((1:Nl)./EEGsim.srate,wintimes);
        
        MODEL = hlp_sift_emptymodel('winStartTimes',winStartTimes ,    ...
                                    'winstep',winStep,'winlen',winLen, ...
                                    'morder',ModelOrder,'algorithm','Oracle');
        if ~all(mu==0)
            MODEL.mu = repmat({mu},1,length(winStartIdx));
        end
        MODEL.PE = repmat({C},1,length(winStartIdx));
        if ~iscell(A)
            MODEL.AR = repmat({A},length(winStartIdx));
        elseif  length(A)==1
            MODEL.AR = repmat(A,length(winStartIdx));
        else
            MODEL.AR = A(winStartIdx+round(winLenPnts/2));
        end
        EEGtrue           = EEGsim;
        EEGtrue.CAT.MODEL = MODEL;
        EEGtrue.setname   = [EEGsim.setname '[Ground Truth]'];
        EEGtrue.condition = [EEGsim.condition '[Ground Truth]'];
        
    end % makeOracle


    % write the simulation code to an m-file
    function writeSim(expr,morder,simName,simDescrip,simAuthor,simRef,verb)

        machineSimName = hlp_variableize(simName);
        
        % generate the function source code
        str = ...
            ['function [expr morder] = sim_ex_' machineSimName '(varargin)\n', ...
            '%% Simulation:  ' simName '\n',            ...
            '%%\n', ...
            '%% Description:  \n', ...  
            '%% \n', ...
            '%% ' simDescrip '\n', ...
            '%% \n', ...
            '%% Recommended Settings: \n', ...
            '%% \n', ...
            '%% Author Credits: \n', ...
            '%% \n', ...
            '%% ' simAuthor '\n', ...
            '%% \n', ...
            '%% References and Code: \n', ...
            '%% \n', ...
            '%% ' simRef '\n', ...
            '%% \n', ...
            '%% ------------------------------------------------------------------------ \n', ...
            '\n', ...
            '%% specify the default system of equations \n', ...
            'expr_def = ' hlp_tostring(expr) ';\n', ...
            '\n', ...
            '%% set up argument definitions \n', ...
            'arg_define(varargin, ... \n', ...
            '    arg({''expr'',''DynamicalEquations''},expr_def,[],''System of equations''), ... \n', ...
            '    arg({''morder'',''ModelOrder''},' num2str(morder) ',[],''Model order. This is mandatory'')); \n', ...
            '\n', ...
            'if isempty(morder) \n', ...
            '    error(''SIFT:sim_examples:badParam'',''ModelOrder must be specified''); \n', ...
            'end \n'];

        % write the code to an m-file
        fname = [fileparts(mfilename('fullpath')) filesep 'examples' filesep 'sim_ex_' machineSimName '.m']; 
        % check if simulation already exists
        if exist(fname,'file')
            res = questdlg2(sprintf('Replace existing simulation %s?',machineSimName),'File Overwrite Warning','Yes','No','No');
            if strcmpi(res,'no')
                return;
            end
        end
        fid   = fopen(fname,'w+');
        if fid==-1
            error('Cannot open file %s for writing',fname); end
        fprintf(fid,str);
        fclose(fid);

        if verb
            fprintf('Simulation saved to %s\n',fname);
        end
        % clear simulation examples cache
        clear hlp_getSimExamples;

    end %writeSim

end % sim_varmodel
    

