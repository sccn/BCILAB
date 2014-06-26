function signal = flt_reconstruct(varargin)
% Reconstruct the given data in a new (possibly overcomplete) basis.
% Signal = flt_reconstruct(Signal, Dictionary, Arguments...)
%
% Sparse Reconstruction is an approach to re-representing a signal in a new basis, which may have
% (far) more basis vectors than the signal has channels [5,7]; in these cases, a unique solution can
% only be found if the additional assumption is incorporated that at any given time there is only a
% small (sparse) set of basis vectors active. Possible use cases include obtaining a sparse
% reconstruction of cortical surface area activation (in terms of a basis of surface patch
% activities).
%
% In:
%   Signal      : epoched or continuous EEGLAB data set
%
%   Dictionary  : [#Channels x #BasisVectors] dictionary of basis vectors in terms of which
%                 the signal shall be reconstructed
%
%   Variant     : reconstruction method to use and parameters; can be one of the following:
%                 'MacKay' : sparse Bayesian Learning using the original MacKay update rule [1,2]
%                 'FastEM' : sparse Bayesian Learning using the fast EM method by Wipf et al. [3]
%                 'TradEM' : sparse Bayesian Learning using a traditional EM method (quite slow) [2]
%                 'FOCUSS' : FOCUSS method for sparse recovery [4]
%                 'l1'     : sparse recovery using standard l1 norm regularization, as in
%                            Compressive Sensing [5]
%                 'tv'     : sparse erecovery using the total-variation norm [6]
%                 'l2'     : simple non-sparse least-squares reconstruction
%                  Note that these variants have optional sub-parameters, which can be passed if 
%                  variant is specified as cell array, e.g., {'MacKay','Lambda',1,'WindowLength',0.5}
%                  (default: 'FastEM')
%
%   BasisInfo   : new channel locations, with one entry per basis vector can either be a cell array 
%                 of channel names or a struct  with field 'labels' and possibly other fields
%                 (default: {'1','2','3',...})
%
%   TransformData : whether to place the result in the .data field of the output signal, or in some
%                   other field (note: currently always true) (default: true)
%
%   Verbose      : whether to show verbose outputs (default: true)
%
% Out:
%   Signal : the signal reconstructed in terms of the dictionary.
%
% Notes:
%   The only parameters that may be specified by position (instead of as name-value pair) are the first two.
%   This function is experimental in nature and currently too slow for online processing.
%
% Examples:
%   % reconstruct the signal in terms of a random dictionary using the default settings
%   eeg = flt_reconstruct(eeg,randn(32,1000))
%
%   % as before, but use the sparse Bayesian Learning using the MacKay update rules
%   eeg = flt_reconstruct(eeg,randn(32,1000),'Variant','MacKay')
%   
%   % as before, but override the default Lambda value to obtain a different sparsity/accuracy tradeoff
%   eeg = flt_reconstruct(eeg,randn(32,1000),'Variant',{'MacKay','Lambda',2})
%
%   % as before, but reconstruct smaller windows at the same time and use more iterations
%   eeg = flt_reconstruct(eeg,randn(32,1000),'Variant',{'MacKay','WindowLength',0.5,'MaxIterations',500})
%
%   % reconstruct the signal in terms of a random dictionary and specify some custom channel labels
%   eeg = flt_reconstruct(eeg,randn(32,1000),'BasisInfo',{'A','B','C','D', ...})
%
%   % reconstruct using l1 norm regularization (standard compressive sensing) and specify a custom noise level
%   eeg = flt_reconstruct(eeg,randn(32,1000),'Variant',{'l1','NoiseLevel',0.05})
%   
%   % use a (non-sparse) least-squares reconstruction
%   eeg = flt_reconstruct(eeg,randn(32,1000),'Variant'l2')
%
% References:
%   [1] MacKay D.J.C. "Bayesian Interpolation"
%       Neural Computation 4(3):415-447 (1992)
%   [2] Tipping M.E. and Smola A., "Sparse Bayesian Learning and the Relevance Vector Machine". 
%       Journal of Machine Learning Research 1: 211?244. (2001)
%   [3] Wipf D.P. and Nagarajan S., "A New View of Automatic Relevance Determination,"
%       In J.C. Platt, D. Koller, Y. Singer, and S. Roweis, editors, Advances in Neural Information Processing Systems 20, MIT Press, 2008.
%   [4] Gorodnitsky I.F., George J.S., Rao B.D., "Neuromagnetic source imaging with FOCUSS: A recursive weighted minimum norm algorithm"
%       J. Electroencephalography and Clinical Neurophysiology, 95(4) (1995).
%   [5] Candes E. J.. "Compressive sampling" 
%       Proceedings of the International Congress of Mathematicians, Madrid, Spain, (2006).
%   [6] Chambolle A., "An algorithm for total-variation minimization and applications"
%       Journal of Mathematical Imaging and Vision, 20 pp. 89-97 (2004)
%   [7] Wipf D.P., Owen J.P., Attias H.T., Sekihara K., and Nagarajan S. "Robust Bayesian Estimation of the Location, Orientation, and Timecourse of Multiple Correlated Neural Sources using MEG" 
%       NeuroImage, vol. 49(1) (2010)
%
% See also:
%   NESTA_UP, sparse_learning
%
% TODO:
%   Upgrade to ICSD method by Ozgur.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-10

% flt_reconstruct_version<0.8> -- for the cache

if ~exp_beginfun('offline') return; end;

% requires relatively clean data, cannot be used as a basis for ICA
declare_properties('name','SparseReconstruction', 'experimental',true, 'follows',{'flt_ica','flt_iir','flt_fir'}, 'independent_channels',false, 'independent_trials',true);

arg_define([0 2],varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'dict','Dictionary'},[],[],'Dictionary of basis vectors. The signal is reconstructed in terms of these basis vectors; formatted as [#Channels x #BasisVectors].','type','expression','shape','row'), ...
    arg_subswitch({'variant','Variant'},'FastEM', {...
        {'MacKay',{...
            arg({'lambda','Lambda'},1,[],'Regularization parameter. Controls the sparsity-accuracy tradeoff.'), ...
            arg({'hyperprior','HyperPrior'},[],[],'Hyper-parameter gamma per basis vector. This is the prior probability of each basis vector. If empty, amounts to the minimum-norm solution.','shape','row'), ...
            arg({'window_len','WindowLength'},1,[],'Window length. Overlapped windows of this length will be decomposed, in seconds.'), ...
            arg({'maxiter','MaxIterations'},100,[],'Maximum number of iterations.')}}, ...
        {'FastEM',{...
            arg({'lambda','Lambda'},1,[],'Regularization parameter. Controls the sparsity-accuracy tradeoff.'), ...
            arg({'hyperprior','HyperPrior'},[],[],'Hyper-parameter gamma per basis vector. This is the prior probability of each basis vector. If empty, amounts to the minimum-norm solution.','shape','row'), ...
            arg({'window_len','WindowLength'},1,[],'Window length. Overlapped windows of this length will be decomposed, in seconds.'), ...
            arg({'maxiter','MaxIterations'},100,[],'Maximum number of iterations.')}}, ...
        {'TradEM',{...
            arg({'lambda','Lambda'},1,[],'Regularization parameter. Controls the sparsity-accuracy tradeoff.'), ...
            arg({'hyperprior','HyperPrior'},[],[],'Hyper-parameter gamma per basis vector. This is the prior probability of each basis vector. If empty, amounts to the minimum-norm solution.','shape','row'), ...
            arg({'window_len','WindowLength'},1,[],'Window length. Overlapped windows of this length will be decomposed, in seconds.'), ...
            arg({'maxiter','MaxIterations'},1000,[],'Maximum number of iterations.')}}, ...
        {'FOCUSS',{...
            arg({'lambda','Lambda'},1,[],'Regularization parameter. Controls the sparsity-accuracy tradeoff.'), ...
            arg({'p','QuasiNorm'},1,[],'p-valued quasi-norm. Smaller values yield sparser solutions One yields the l1 norm.'), ...
            arg({'window_len','WindowLength'},1,[],'Window length. Overlapped windows of this length will be decomposed, in seconds.'), ...
            arg({'maxiter','MaxIterations'},10,[],'Maximum number of iterations.')}}, ...
        {'l1',{...
            arg({'noiselev','NoiseLevel','sigma'},0.1,[],'Noise std dev estimate.'), ...
            arg({'mufinal','Tolerance'},[],[],'Solution accuracy tolerance. If empty, a heuristic will be used.','shape','scalar'), ...
            arg({'maxiter','MaxIterations'},3000,[],'Maximum number of iterations.'), ...
            arg({'maxintiter','MaxContinuations'},5,[],'Maximum number of continuation steps.')}}, ...
        {'l2',{}}}, 'Optimization method to use. MacKay is Sparse Bayesian Learning with MacKay update rules (from the paper Bayesian Interpolation), FastEM is a fast Expectation-Maximization (EM) based update rule (from David Wipf), TradEM is the traditional EM rule (from Mike Tipping''s SBL papers), FOCUSS is the FOCUSS algorithm from Bhaskar Rao, Kreutz-Delgado and Gorodnitsky, l1 is the traditional l1-norm recovery, tv is recovery using the total variation norm, and l2 is a non-sparse pseudoinverse-based approach.'), ...            
    arg({'basisinfo','BasisInfo'},[],[],'Chanlocs for new signal basis. This is either a cell array of channel names or a struct array with a field ''labels''. If unspecified, defaults to {''1'',''2'',''3'',...}.','type','expression','shape','row'), ...
    arg({'dotransform','TransformData'},true,[],'Transform data into new representation. If false, the data set will be annotated with the field sparseact, which contains the activity.'), ...
    arg({'verbose','Verbose'},true,[],'Report progress of algorithm.','cat','Miscellaneous'));

[C,S,T] = size(signal.data);
if isempty(dict) 
    dict = eye(C); end

% reconstruct data
signal.data = reshape(signal.data,C,[]);
switch lower(variant)
    case 'l1'
        % use l1 reconstruction
        data = double(signal.data);
        signal.data = zeros(size(dict,2),S*T);
        try
            % use NESTA
            dictnorm = norm(dict*dict');            
            lambda = variant.noiselev * sqrt(2*log(C));
            if isempty(variant.mufinal)
                variant.mufinal = 0.1*variant.noiselev/dictnorm; end            
            for t=1:S*T
                signal.data(:,t) = NESTA_UP(dict,[],data(:,t),lambda,dictnorm,variant.mufinal,struct('MaxIntIter',variant.maxintiter,'Verbose',verbose,'TypeMin',hlp_rewrite(variant,'l1','L1'))); end
        catch
            % use the CVX fallback
            N = size(dict,2);
            for t=1:S*T
                cvx_begin
                    variables X(N)
                    minimize(norm(X,1))
                    subject to
                        norm(dict*X - data(:,t)) <= variant.noiselev;
                cvx_end
                signal.data(:,t) = X;
            end
        end
    case {'mackay','fastem','tradem','focuss'}        
        % use sparse Bayesian learning
        mode = struct('mackay',0,'fastem',1,'tradem',2,'focuss',[3 variant.p]);
        if isempty(variant.hyperprior)
            variant.hyperprior = 0; end        
        % process the signal in overlapped windows
        window_len = variant.window_len*signal.srate;
        wnd = 0:window_len-1;
        wnd_weight = repmat(hann(length(wnd))',C,1);
        offsets = 1 + floor(0:window_len/2:(S*T)-window_len);
        W = length(offsets);
        X = signal.data;
        signal.data = zeros(size(dict,2),S*T);
        for o=1:W
            S = X(:,offsets(o) + wnd).*wnd_weight;
            signal.data(:,offsets(o)+wnd) = signal.data(:,offsets(o)+wnd) + sparse_learning_fast(dict,S,variant.lambda,variant.maxiter,mode.(variant),variant.hyperprior,variant.verbose);
        end
    case 'l2'
        % use the pseudoinverse
        signal.data = dict\signal.data;
    otherwise
        error(['unsupported variant selected: ' variant]);
end

if dotransform
    % override data
    signal.data = reshape(signal.data,[],S,T);
    signal.nbchan = size(signal.data,1);
    % rewrite chanlocs
    if isempty(basisinfo)
        signal.chanlocs = struct('labels',cellfun(@num2str,num2cell(1:signal.nbchan,1),'UniformOutput',false));
    elseif length(basisinfo) == signal.nbchan
        if isfield(basisinfo,'labels')
            signal.chanlocs = basisinfo;
        elseif iscellstr(basisinfo)
            signal.chanlocs = struct('labels',basisinfo);
        else
            error('unsupported format for the ''basisinfo'' parameter');
        end        
    else
        error('length of ''basisinfo'' parameter does not match number of basis vectors');
    end
else
    error('dotransform is currently required to be 1');
end

exp_endfun;
