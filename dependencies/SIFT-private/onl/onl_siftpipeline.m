function [signal g] = onl_siftpipeline(varargin)
% implements a basic SIFT pipeline (Preprocessing, Modeling, Connectivity)

% if ~exp_beginfun('filter'), return; end
% declare_properties('name','SIFT');

try
    EEG = arg_extract(varargin,{'EEG','signal','Signal'},struct([]));
catch e
    EEG = struct([]);
end

g = arg_define(varargin, ...
        arg_norep({'EEG','Signal','signal'}), ...
        arg({'channels','Channels'}, [], [], 'Cell array of channel names to retain.','type','cellstr','shape','row'), ...
        arg_sub({'preproc','Preprocessing'},{'EEG',EEG},@pre_prepData,'Pre-processing options'), ...
        arg_subswitch({'modeling','Modeling'},{'Segmentation VAR' 'EEG',EEG}, ...
            hlp_getModelingApproaches, 'Select a modeling approach and define parameters.'), ...
        arg_subtoggle({'selModelOrder','AutoSelectModelOrder'},{},...
        { ...  % FIXME: for some reason setting this as defaults in modelOrderSelection gens an error: {'plot','off','icselector',{'hq'}}
            arg_sub({'modelOrderSelection','ModelOrderSelection'},{'plot','off','icselector',{'hq'}},@est_selModelOrder,'Model order selection control.','suppress',{'modelingApproach'}), ...
            arg_subswitch({'minimizer','Selector'},{'min'}, ...
            {'elbow', {}, ...
             'min',   {}, ...
             'prctile' ...
             {arg({'prclim','PercentileLimits'},90,[1 100],'Upper percentile limit for order selection. If PercentileLimits = L, This selects the model order, p, for which L% of all sampled windows indicate an optimal model order of p or lower.')} ...
             },'Method for determining optimal model order. If "min", then the optimal model order is the one that minimizes the information criterion. If "elbow" then the optimal order is the elbow of the function of informaton criterion versus model order. If "prctile" then percentile is used.'), ...
        },'Automatically select model order'), ...
        arg_subtoggle({'connectivity','Connectivity'},{'EEG',EEG,'MODEL',struct([])},@est_mvarConnectivity,'Select connectivity methods'), ...
        arg_subtoggle({'validation','Validation'},{'EEG',EEG},@est_validateMVAR,'Validate Model fit'), ...
        arg({'printValidation','PrintValidation'},false,[],'Print validation output to console'), ...
        arg({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
        );
clear EEG;

if g.selModelOrder.arg_selection && length(g.selModelOrder.modelOrderSelection.icselector) > 1
    error('SIFT:onl_siftpipeline:MoreThanOneIC','Only one InformationCriterion can be selected for ModelOrderSelection.'); end
            

% Implement basic SIFT pipeline (Preprocessing, Modeling, Connectivity)
% -------------------------------------------------------------------------

if ~isempty(g.channels) && ~isequal(g.channels,1:g.EEG.nbchan) && ~isequal(g.channels,{char(zeros(1,0))}) % TODO: remove last isequal hack
    if g.verb
        g.EEG = pop_select(g.EEG,'channel',g.channels,'sorttrial','off'); 
    else
        [console_text,g.EEG] = evalc('pop_select(g.EEG,''channel'',g.channels,''sorttrial'',''off'');');  %#ok<ASGLU>
    end
end

% rereference to average
% g.EEG.data = bsxfun(@minus,g.EEG.data,mean(g.EEG.data)); end

% pre-process data
g.EEG = pre_prepData('EEG',g.EEG,g.preproc,'verb',g.verb,'arg_direct',true);

% if strcmp(g.preproc.sigtype.arg_selection,'Sources')
%     % if using sources, construct the dipfit matrix
%     g.EEG.dipfit = hlp_microcache('dipfit',@hlp_makeDipfitStruct,hmObj.sourceSpace,g.EEG.roiVertices);
% end
                
% get the m-file name of the function implementing the modeling approach
modelingFuncName = hlp_getModelingApproaches('mfileNameOnly',g.modeling.arg_selection);

% determine optimal model order
if g.selModelOrder.arg_selection
    % model order selection
    
    g.selModelOrder.modelOrderSelection.modelingApproach = g.modeling;

    IC = est_selModelOrder('EEG',g.EEG,g.selModelOrder.modelOrderSelection,'arg_direct',true);
    if iscell(IC), IC = IC{1}; end
    
    icselector = g.selModelOrder.modelOrderSelection.icselector{1};
    switch g.selModelOrder.minimizer.arg_selection
        case 'min'
            popt = ceil(mean(IC.(icselector).popt));                    %  mean
        case 'elbow'
            popt = ceil(mean(IC.(icselector).pelbow));
        case 'prctile'
            popt = ceil(prctile(IC.(icselector).popt,g.selModelOrder.minimizer.prclim));   %  upper 95th prctile
    end
    
    g.modeling.morder = popt;
end
            
% fit model
g.EEG.CAT.MODEL = feval(modelingFuncName,'EEG',g.EEG,g.modeling,'verb',g.verb,'arg_direct',true);

% perform model validation (optional)
if g.validation.arg_selection
    [whitestats PCstats stability residualstats] = est_validateMVAR('EEG',g.EEG,g.validation,'arg_direct',true);

    if isempty(whitestats)
        whitenessCriteria = {};
    else
        whitenessCriteria = g.validation.checkWhiteness.whitenessCriteria;
    end
    
    if g.printValidation
        vis_validation(whitestats,PCstats,stability,whitenessCriteria);
    end
    
    % store results
    g.EEG.CAT.validation.whitestats = whitestats;
    g.EEG.CAT.validation.PCstats    = PCstats;
    g.EEG.CAT.validation.stability  = stability;
    g.EEG.CAT.validation.residualstats = residualstats;
end

% calculate connectivity
if g.connectivity.arg_selection
    g.EEG.CAT.Conn = est_mvarConnectivity('EEG',g.EEG,'MODEL',g.EEG.CAT.MODEL,g.connectivity,'verb',g.verb,'arg_direct',true);
else
    g.EEG.CAT.Conn = [];
end

% return SIFT-augmented dataset
signal = g.EEG;

if nargout>1
    g = rmfield(g,'EEG');
end