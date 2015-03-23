function [W,S,M,initW] = rica(X,params,options)
% Perform Independent Component Analysis via Reconstruction ICA (RICA)
% [Weights,Sphere,ChanList,InitWeights] = rica(X,Parameters,SolverOptions)
%
% In:
%   X : the data [#channels x #samples]
%
%   Parameters : structure of options with fields:
%       .lambda : regularization parameter that controls the strength of the
%                 sparsity / independence cost (default: 0.05)
%       .gamma  : regularization parameter that controls the strength of the
%                 anatomical constraints, if any (default: 0.01)
%       .numFeatures : number of components to learn (default: #channels)
%                      this can be larger than the number of channels
%       .epsilon : smoothing parameter for sparsity cost near minimum (default: 1e-5)
%       .dict_criterion : criterion for dictionary learning (can be:
%                         * 'reconstruction' for vanilla RICA,
%                         * 'cortically_anchored' for RICA with constraint term);
%                           (default: 'reconstruction')
%       .initialization : how to initialize the components (can be:
%                         * 'radial' for radially oriented components
%                         * 'random' for randomly initialized components
%       .cov_blocksize : blocksize for calculating the covariance robustly;
%                        robustness to shorter-than-blocksize samples artifacts cannot be
%                        guaranteed (but uses less memory for larger sizes);
%                        (default: 10)
%       .random_init_scale : level of randomness to initial value (default: 0.1)
%       .cov_rejection : reject data samples beyond this number of standard deviations
%                        from the robustly estimated data distribution (default: 5)
%       .temporal_normalization : whether to normalize each sample (default: false)
%       .anchors : optionally a cell array of anchor subspaces;
%                  each is either a matrix of #channels x #dimension forward projections
%                  (these are typically x/y/z orthogonal dipolar projections emitted
%                   from a given location, where each cell contains a different location)
%                  or a vector of anatomical labels (see LF_centroidmap in code below)
%       .anchor_init : how anchors shall be initialized (can be:
%                      * 'from_subspace' : anchors are initialized to randomized vectors in the
%                                          respective subspace
%                      * 'bestmatch' : anchored components are initialized with the best-matching
%                                      components out of the initial weights (best match = lowest constraint violation)
%                      * 'perpendicular' : perpendicular to cortex
%                      * cell array of initial vectors
%       .randseed : random seed of initial solution (default: 10)
%       .max_restarts : maximum number of restarts (if weights blow up) (default: 20)
%       .topoplot : whether to plot progress in form of topoplots (default: false)
%       .chan_labels : list of 10-20 channel labels (for anatomical lookup support)
%
%  SolverOptions: structure of options with fields (see also minFunc):
%       .Method : optimization method to use (default: 'scg')
%       .MaxIter : maximum number of iterations (default: 150)
%       .Display : what to display in the console; can be 'none','iter', or 'full' (default: 'iter')
%
% Out:
%   Weights : unmixing matrix in sphered space
%   Sphere : sphering matrix
%   ChanList : list of retained channel indices


% configure cost function
if ~isfield(params,'lambda')
    params.lambda = 0.05; end
if ~isfield(params,'gamma')
    params.gamma = 0.01; end
if ~isfield(params,'theta')
    params.theta = 1; end
if ~isfield(params,'epsilon')
    params.epsilon = 1e-5; end
if ~isfield(params,'dict_criterion')
    params.dict_criterion = 'reconstruction'; end % 'reconstruction' or 'cortically_anchored'
if ~isfield(params,'initialization')
    params.initialization = 'radial'; end % 'radial', 'random', or cell array of initial vectors
if ~isfield(params,'cov_blocksize')
    params.cov_blocksize = 10; end
if ~isfield(params,'random_init_scale')
    params.random_init_scale = 0.1; end
if ~isfield(params,'cov_rejection')
    params.cov_rejection = 5; end
if ~isfield(params,'temporal_normalization')
    params.temporal_normalization = false; end
if ~isfield(params,'chan_labels')
    params.chan_labels = []; end
if ~isfield(params,'anchors')
    params.anchors = {}; end
if ~isfield(params,'anchor_init')
    params.anchor_init = 'from_subspace'; end % 'from_subspace', 'bestmatch' (={}), 'perpendicular', or cell array of vectors
if ~isfield(params,'randseed')
    params.randseed = 10; end
if ~isfield(params,'max_restarts')
    params.max_restarts = 20; end
if ~isfield(params,'topoplot')
    params.topoplot = false; end
if ~isfield(params,'numFeatures')
    params.numFeatures = []; end

% configure basic optimization parameters
if ~isfield(options,'Method')
    options.Method = 'scg'; end
if ~isfield(options,'MaxIter')
    options.MaxIter = 150; end
if ~isfield(options,'Display')
    options.Display = 'iter'; end

% support for string-formatted anatomical regions (requires 10-20 labels and Guido Nolte's forward modeling toolbox)
if ~isempty(params.anchors) && iscellstr(params.anchors)
    [params.chan_labels,params.anchors,params.anchor_init,channel_mask] = hlp_diskcache('headmodel',@calc_anatomical_constraints,params.chan_labels,params.anchors,params.anchor_init);
    X = X(channel_mask,:);
    M = channel_mask;
end

params.n = size(X,1);
if isempty(params.numFeatures)
    params.numFeatures = size(X,1); end

%% calc robust Gaussian stats
Mu = median(X,2);
Sigma = hlp_diskcache('filterdesign',@cov_blockgeom,X',params.cov_blocksize);
params.Sigma = Sigma;
%% Optionally reject bad samples
if params.cov_rejection && exist('reject_samples_robcov','file')
    X = reject_samples_robcov(X,params.cov_rejection,Mu,Sigma); end
%% Whiten the data
X = bsxfun(@minus,X,Mu);
mix = sqrtm(Sigma);
params.mix = mix;
sphere = inv(mix);
X = sphere*X; %#ok<*MINV>

%% optionally normalize temporally
if params.temporal_normalization
    m = sqrt(sum(X.^2) + (1e-8));
    X = bsxfun(@rdivide,X,m);
end

params.sphere = sphere;
params.mix = mix;
params.rcov = Sigma;

%% Select cost function
if strcmp(params.dict_criterion,'reconstruction')
    optFunc = @(theta) softICACost(theta, X, params);
elseif strcmp(params.dict_criterion,'cortically_anchored')
    I = eye(size(X,1));
    % pre-calculate the subspace matrices...
    if isempty(params.anchors)
        error('You need to specify a .anchors field that is a cell array of Nxd (usually d=3) dipolar forward projection axes, one cell per location.'); end
    for l = 1:length(params.anchors)
        U = sphere*params.anchors{l};
        % renormalize
        U = bsxfun(@times,U,1./sqrt(sum(U.^2)));
        if size(U,1) ~= length(I)
            error('The anchors must be [#channels x #subspace-dimension].'); end
        % create the projection weighting matrix and bake in the regularization parameter, too
        params.subspaces{l} = (I - U*inv(U'*U)*U') * params.gamma; %#ok<MINV>
    end
    optFunc = @(theta) softICACost3(theta, X, params);
elseif strcmp(params.dict_criterion,'incoherence')
    optFunc = @(theta) softICACost2(theta, X, params);
elseif strcmp(params.dict_criterion,'simple')
    optFunc = @(theta) softICACost_simple(theta, X, params);
else
    error('Unknown optimization criterion!');
end

%% set output function if desired
if iscell(params.topoplot) || params.topoplot
    if ~iscell(params.topoplot)
        params.topoplot = {}; end
    f = figure(params.topoplot{:});
    inputhash = hlp_fingerprint({X,params});
    options.outputFcn = @(x,varargin)display_output(x,params,f,inputhash,varargin{:});
end

%% Init weights
randn('state',params.randseed); %#ok<RAND>
initTheta = init_weights(params,sphere);
%opttheta = initTheta;

%% Run optimization; restart if weights blow up
for k=1:params.max_restarts
    try
        [opttheta,cost,exitflag] = minFunc(optFunc,initTheta(:),options); %#ok<ASGLU>
        if exitflag<0
            error('rica:abnormalTermination','Abnormal termination.');
        else
            break;
        end
    catch e
        if strcmp(e.identifier,'rica:abnormalTermination')
            disp('Weights blew up or got stuck; restarting with new randomized weights.');
            initTheta = init_weights(params,sphere);
            continue;
        else
            rethrow(e);
        end
    end
end

%% finalize results
W = reshape(opttheta,params.numFeatures,params.n);
S = sphere;
if ~exist('M','var')
    M = true(1,params.n); end
M = find(M);
initW = initTheta;
end

% debug display:
% global debug_chanlocs;figure;topoplot_grid(pinv(W*S),debug_chanlocs)


function initTheta = init_weights(params,sphere)
% initialize randomly
if strcmp(params.initialization,'random_dipoles')
    atlas = io_load('resources:/sa_montreal_small.mat');
    % get a matrix of random dipoles...
    initTheta = [];
    for k=params.numFeatures:-1:1
        voxel = 1+floor(rand*(size(atlas.sa.V_medium,2)-1));
        orientationWeights = randn(3,1);
        orientationWeights = orientationWeights/norm(orientationWeights);
        initTheta(k,:) = squeeze(atlas.sa.V_medium(:,voxel,:))*orientationWeights;
    end
    % get the chanlocs of the head model (10-20 system)
    refLocs = getfield(io_load('resources:/caps/Standard-10-5-Cap385.cap','-mat'),'CAP');
    [dummy,atlasLocIdx,atlasSubset] = intersect(lower({refLocs.labels}),lower(atlas.sa.clab_electrodes)); %#ok<ASGLU>
    oldLocs = [refLocs(atlasLocIdx).X; refLocs(atlasLocIdx).Y; refLocs(atlasLocIdx).Z];
    % get the chanlocs of the data 
    if any([length([params.chan_locs.X]),length([params.chan_locs.Y]),length([params.chan_locs.Z])] < length(params.chan_locs))
        error('Some of your channels have no locations; cannot use the random_dipoles initialization of rica.'); end
    % interpolate initial weights from head model locs to data locs
    newLocs = [params.chan_locs.X; params.chan_locs.Y; params.chan_locs.Z];
    interpMat = sphericalSplineInterpolate(oldLocs,newLocs);
    initTheta = initTheta(:,atlasSubset)*interpMat';
    % handle sphering
    initTheta = initTheta/params.Sigma;
    initTheta = initTheta*inv(sphere);
else
    initTheta = strcmp(params.initialization,'radial') * eye(params.numFeatures,params.n) + randn(params.numFeatures,params.n)*(params.random_init_scale+strcmp(params.initialization,'random'));
end
initTheta = initTheta ./ repmat(sqrt(sum(initTheta.^2,2)), 1, size(initTheta,2));

if strcmp(params.dict_criterion,'cortically_anchored')
    if iscell(params.anchor_init) && ~isempty(params.anchor_init)
        if length(params.anchor_init) == params.numFeatures
            initTheta = [];
            for l=1:length(params.subspaces)
                initTheta(l,:) = sphere*params.anchor_init{l}; end
        else
            % partially override the initialization if desired
            for c=1:length(params.anchor_init)
                % sphere and renormalize the anchor component
                comp = sphere*params.anchor_init{c};
                comp = comp/norm(comp);
                % replace the formerly best-matching component from current initialization
                match = abs(sum(bsxfun(@times,initTheta',comp)));
                [dummy,idx] = max(match); %#ok<ASGLU>
                % switch with the c'th component and override the c'th component
                initTheta(idx,:) = initTheta(c,:);
                initTheta(c,:) = comp;
            end
        end
    elseif strcmp(params.anchor_init,'bestmatch') || isempty(params.anchor_init)
        % reorder initial components based on lowest cost under the constraints
        for l=1:length(params.subspaces)
            for n=1:params.numFeatures
                cost(n) = initTheta(n,:) * params.subspaces{l} * initTheta(n,:)'; end
            [dummy,bestidx] = min(cost); %#ok<ASGLU>
            % put the component there
            tmp = initTheta(l,:);
            initTheta(l,:) = initTheta(bestidx,:);
            initTheta(bestidx,:) = tmp;
        end
    elseif strcmp(params.anchor_init,'from_subspace')
        if length(params.subspaces) == params.numFeatures
            % use the components as they are
            initTheta = [];
            for l=1:length(params.subspaces)
                initTheta(l,:) = sphere*params.anchors{l}*randn(size(params.anchors{l},2),1); end
        else
            for l=1:length(params.subspaces)
                sp = sphere*params.anchors{l};
                comp = sp*randn(size(sp,2),1);
                comp = comp/norm(comp);
                initTheta(l,:) = comp;
            end
        end
    else
        error('Unknown anchor initialization');
    end
end

% renormalize again, for good measure
initTheta = initTheta ./ repmat(sqrt(sum(initTheta.^2,2)), 1, size(initTheta,2));
%global debug_chanlocs;figure;topoplot_grid(inv(sphere)^2*(initTheta'*sphere),debug_chanlocs)
%global debug_chanlocs;figure;topoplot_grid((initTheta(1:63,:)*sphere)*inv(sphere)^2,debug_chanlocs)
end

function stop=display_output(x,params,fig,inputhash,mode,iter,varargin)
clf;
x=reshape(x,params.numFeatures,params.n);
topoplot_grid((x(1:end-1,:)*params.sphere)*inv(params.sphere)^2,hlp_microcache('chanlabels',@set_infer_chanlocs,params.chan_labels));
drawnow;
save2pdf(sprintf('rica_%i_frame_%i.png',inputhash,iter),gcf);
stop = false;
end

function [LF_chans,LF_anchors,LF_anchorinit,channel_mask] = calc_anatomical_constraints(LF_chans,LF_anchors,LF_anchorinit)
% lookup table of anatomical region centroids (generated by MoBIlab)
% warning: these are not necessarily coregistered particularly well with the head model
LF_centroidmap = struct('Precentral_L', {[4.23108288452151 0.305597157992957 5.12657838054604]}, ...
    'Precentral_R', {[-4.20528020132841 0.673804036332546 5.12761812250973]}, 'Frontal_Sup_L', ...
    {[2.08858531558692 -3.88674676095795 4.19641636320092]}, 'Frontal_Sup_R', {[-2.21802251856876 ...
    -3.42197526707012 4.23797390857122]}, 'Frontal_Sup_Orb_L', {[1.670197560078 -4.95343367704499 ...
    -1.91436767127854]}, 'Frontal_Sup_Orb_R', {[-1.5439059147638 -4.6981967612907 -2.1290543308561]}, ...
    'Frontal_Mid_L', {[3.44601607706013 -3.40935856243948 3.58258171502308]}, 'Frontal_Mid_R', ...
    {[-3.65848090239564 -3.31735265619251 3.44006208905955]}, 'Frontal_Mid_Orb_L', {[3.27423289629387 ...
    -5.53036266593331 -1.24890075352431]}, 'Frontal_Mid_Orb_R', {[-3.15007344473511 -5.58049284721048 ...
    -1.37598923793298]}, 'Frontal_Inf_Oper_L', {[5.57607765119148 -1.22160092601595 1.87820668235864]}, ...
    'Frontal_Inf_Oper_R', {[-5.12518654546578 -1.51103974148078 2.2835043538629]}, 'Frontal_Inf_Tri_L', ...
    {[5.25720954524236 -2.94214626022688 1.4181372407483]}, 'Frontal_Inf_Tri_R', {[-4.99288004332893 ...
    -2.99128831940483 1.37548200375666]}, 'Frontal_Inf_Orb_L', {[4.30081007236675 -3.41241828575515 ...
    -1.44778170896923]}, 'Frontal_Inf_Orb_R', {[-4.39691838955149 -3.32260997520098 ...
    -1.48797082820651]}, 'Rolandic_Oper_L', {[5.71451252012493 -0.111240598724779 0.557889853517538]}, ...
    'Rolandic_Oper_R', {[-5.25432022222563 0.327809974660407 0.927191197017338]}, 'Supp_Motor_Area_L', ...
    {[0.431698144489469 -0.587840261923926 6.06598813456952]}, 'Supp_Motor_Area_R', ...
    {[-0.552551234546444 -0.344963213444904 6.43604620799783]}, 'Frontal_Sup_Medial_L',{[0.339386309547571 -5.03922302451318 3.05529314773942]}, 'Frontal_Sup_Medial_R', ...
    {[-0.599042403013505 -4.95358548975335 3.3315470608752]}, 'Frontal_Med_Orb_L', {[0.553159296736419 ...
    -5.28242672178415 -0.860062252760082]}, 'Frontal_Med_Orb_R', {[-0.52237222952882 -5.28843478796661 ...
    -0.847506190494953]}, 'Insula_L', {[4.04000316273526 -0.991626367382352 0.182033060791309]}, ...
    'Insula_R', {[-4.20804313472544 -0.599463076265843 0.117060827883777]}, 'Cingulum_Ant_L', ...
    {[0.158998122628773 -3.04191418890571 1.60976310192608]}, 'Cingulum_Ant_R', {[-0.341672612264182 ...
    -3.51016844089611 1.36841411085656]}, 'Cingulum_Mid_L', {[0.323654032129958 1.2701467677249 ...
    4.09982437558004]}, 'Cingulum_Mid_R', {[-0.339759220878602 1.22666761520706 3.80942028260427]}, ...
    'Cingulum_Post_L', {[0.253799762446454 3.80270573043469 2.18616811870665]}, 'Cingulum_Post_R', ...
    {[-0.295905097743932 3.92068635848129 1.80072867801197]}, 'Hippocampus_L', {[2.32787239524304 ...
    1.8149338343237 -1.49492324132927]}, 'Hippocampus_R', {[-2.73086294085504 1.5871874600083 ...
    -1.62667625879869]}, 'ParaHippocampal_L', {[1.61074664502006 1.12367086624508 -2.21356921544053]}, ...
    'ParaHippocampal_R', {[-1.70007755539167 1.0863331955891 -2.2987395682162]}, 'Calcarine_L', ...
    {[0.644743908779623 7.86638680536702 0.41530545769496]}, 'Calcarine_R', {[-1.34655674522778 ...
    7.0622525121614 0.815612927283389]}, 'Cuneus_L', {[0.436174227506558 8.0521105184353 ...
    2.75411253307904]}, 'Cuneus_R', {[-1.13632632141353 7.79039576946586 2.9509700661813]}, ...
    'Lingual_L', {[1.45764022312239 7.19528124396283 -0.983516476949735]}, 'Lingual_R', ...
    {[-1.53027886287718 6.8088926535018 -0.80693847227352]}, 'Occipital_Sup_L', {[1.727720655485 ...
    8.91038621986873 3.1231692916246]}, 'Occipital_Sup_R', {[-2.48782529805265 8.60494144930644 ...
    3.01411154849547]}, 'Occipital_Mid_L', {[3.58539669462852 8.47122369159517 1.54439747300616]}, ...
    'Occipital_Mid_R', {[-3.82169076221469 8.403165188002 1.87962910627168]}, 'Occipital_Inf_L', ...
    {[4.32560067799949 8.05473703169455 -1.03584254747973]}, 'Occipital_Inf_R', {[-4.13915120052877 ...
    8.66732901750285 -0.972759221390047]}, 'Fusiform_L', {[3.03306995106131 3.88244101945737 ...
    -2.63904244654184]}, 'Fusiform_R', {[-3.09453956058322 3.62827635080648 -2.58603675436737]}, ...
    'Postcentral_L', {[4.63010409497207 2.01440498124974 4.80059572115281]}, 'Postcentral_R', ...
    {[-4.3993780915274 2.43209019429645 5.0528328634925]}, 'Parietal_Sup_L', {[2.56479927818908 ...
    6.03806817199842 6.06260183837795]}, 'Parietal_Sup_R', {[-2.61847877039288 6.02781590621746 ...
    6.12049304611679]}, 'Parietal_Inf_L', {[4.51081835603749 4.46688051185214 4.70072181090452]}, ...
    'Parietal_Inf_R', {[-4.54377179301509 4.56986897251545 4.96200488139722]}, 'SupraMarginal_L', ...
    {[5.94198387433377 3.5205327966889 3.16734349409244]}, 'SupraMarginal_R', {[-6.04904888223408 ...
    3.29900208232786 3.36174439727687]}, 'Angular_L', {[4.80369100364743 6.24353421545146 ...
    3.54538534104306]}, 'Angular_R', {[-4.77806501884985 6.11227265248037 3.75232035290023]}, ...
    'Precuneus_L', {[0.563440286763802 5.6826217224422 5.02317758707648]}, 'Precuneus_R', ...
    {[-0.44343372640072 5.62909663459633 4.82196678692494]}, 'Paracentral_Lobule_L', ...
    {[0.540521246467874 2.51820695734066 6.73651427005108]}, 'Paracentral_Lobule_R', ...
    {[-0.309626690701906 3.30773124242547 7.15953489872043]}, 'Temporal_Sup_L', {[5.53356669307572 ...
    1.70437586710982 0.470132860954377]}, 'Temporal_Sup_R', {[-5.88637424887785 2.06929301884502 ...
    0.846201776642597]}, 'Temporal_Pole_Sup_L', {[4.01021113499556 -1.67284654220415 ...
    -2.14466608308187]}, 'Temporal_Pole_Sup_R', {[-4.44289359026018 -1.64435707946637 ...
    -1.85674374989261]}, 'Temporal_Mid_L', {[6.03983234828478 3.621427188918 -0.108595187168253]}, ...
    'Temporal_Mid_R', {[-6.06270905130349 3.8417450483484 -0.0829066147782176]}, 'Temporal_Pole_Mid_L', ...
    {[3.8384686505636 -1.58861587930798 -3.89716162887448]}, 'Temporal_Pole_Mid_R', {[-4.68551379712917 ...
    -1.709265911014 -3.34739979850841]}, 'Temporal_Inf_L', {[5.46129029296911 2.85901736886324 ...
    -2.87284412264196]}, 'Temporal_Inf_R', {[-5.53957465355057 2.91863205101178 -2.77216211203766]}, ...
    'Olfactory_L', {[0.17374500963458 -1.07905767596242 -1.18429284451753]}, 'Olfactory_R', ...
    {[-0.44084981073756 -1.41739466639547 -0.969831087606087]}, 'Rectus_L', {[0.29754183833035 ...
    -3.95988559434341 -2.28624398338371]}, 'Rectus_R', {[-0.511346423329275 -3.57672239270401 ...
    -2.33852965006135]}, 'Amygdala_L', {[2.46903665702496 -0.339323993784056 -1.62164758253021]}, ...
    'Amygdala_R', {[-2.46903665702496 -0.339323993784056 -1.62164758253021]}, 'Caudate_L', ...
    {[0.0323034230426779 -0.959211470373974 0.625028875228461]}, 'Caudate_R', {[-0.0787586435037352 ...
    -1.48085818901111 0.146815044083548]}, 'Thalamus_L', {[0.20276678183809 1.4463228836034 ...
    0.661798544960814]}, 'Thalamus_R', {[-0.260684002425553 1.53329231645962 0.342886396190338]}, ...
    'Heschl_L', {[4.56872759776798 1.61039340422665 0.764584372615752]}, 'Heschl_R', ...
    {[-4.72127249282475 1.59334403525277 1.01890557792087]});


if isempty(LF_chans)
    error('You can pass in anatomical ROI anchor labels only if you also specify 10-20 channel labels.'); end
fprintf('Calculating anatomical constraint parameters from head model...');
% initialize source analysis toolbox
sa = prepare_sourceanalysis(LF_chans, 'montreal');
% reduce missing channels
[dummy,missing_channels] = setdiff(lower(LF_chans),lower(sa.clab_electrodes)); %#ok<ASGLU>
channel_mask = true(1,length(LF_chans));
channel_mask(missing_channels) = false;
LF_chans = LF_chans(channel_mask);
% incorporate common average reference (but note that the ICA might run into trouble with that...)
sa.fp.lintrafo = eye(length(LF_chans)) - ones(length(LF_chans))/length(LF_chans);
% find closest vertices in fine (in-cortex) mesh
for c=1:length(LF_anchors)
    LF_centroids(c,:) = LF_centroidmap.(LF_anchors{c});
    [mindiff,bestidx] = min(sqrt(sum(bsxfun(@minus,sa.grid_fine_incortex,LF_centroids(c,:)).^2,2))); %#ok<ASGLU>
    LF_centroids_grid(c,:) = sa.grid_fine_incortex(bestidx,:);
end

% also look up the normal for that location
for c=1:size(LF_centroids,1)
    [mindiff,bestidx] = min(sqrt(sum(bsxfun(@minus,sa.cortex.vc,LF_centroids(c,:)).^2,2))); %#ok<ASGLU>
    LF_normals(c,:) = sa.normals_cortex(bestidx,:);
end

% generate forward projections
for c=1:size(LF_centroids_grid,1)
    vectors = {[1 0 0],[0 1 0],[0 0 1]};
    for ax = 1:length(vectors)
        LF(:,c,ax) = forward_general([LF_centroids_grid(c,:),vectors{ax}], sa.fp); end
    LF_normfield(:,c) = forward_general([LF_centroids_grid(c,:),LF_normals(c,:)], sa.fp);
end

% override anchor parameters
for c=1:size(LF,2)
    LF_anchors{c} = squeeze(LF(:,c,:)); end
if strcmp(LF_anchorinit,'perpendicular')
    LF_anchorinit = [];
    for c=1:size(LF_normfield,2)
        LF_anchorinit{c} = squeeze(LF_normfield(:,c)); end
end
disp('done.');
end