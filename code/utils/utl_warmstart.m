function model = utl_warmstart(initial_model,learning_function,parameter_index,varargin)
% Utility function to ease the definition of warm-start techniques in BCILAB.
% Model = utl_warmstart(InitialModel,LearningFunction,SpecialParamIndex,Arguments...)
%
% Warm-starting is an approach in which a result is computed based on some previously computed result,
% potentially saving redundant computations. 
% 
% Warm-starting is frequently used to solve a sequence of constrained optimization problems with
% incrementally larger (i.e. less strict) constraint sets, where the solution of the previous run
% is used as initial value for the next run. In certain regularized machine learning settings, 
% it can be used to efficiently learn the entire regularization path for a series of regularization 
% parameter values.
%
% This function simplifies this task in BCILAB by automatically keeping track of the model that was 
% previously computed for a given parameter setting (but possibly for another regularization 
% parameter value), and using it (in place of the InitialModel) as parameter to the LearningFunction.
% Note that when implementing the parameter search loop, cross-validation, etc. ad hoc by hand, 
% there is little need for this function - it only helps when it is difficult or impossible to 
% manually keep track of the previous models.
%
% To be able to distinguish between the learner arguments (out of Arguments) that vary over a
% sequence of related learning problems (namely the regularization parameters) and those that
% constitute a different (unrelated) learning problem (such as the training data itself), the
% function requires that the index (or indices, into Arguments) of those parameters that may vary is
% passed as SpecialParamIndex.
%
% In:
%   InitialModel : the initial model to use for the learning function
%
%   LearningFunction : function handle to learn a new model given a previous model and some 
%                      arbitrary parameters; invoked as NewModel = LearningFunction(PrevModel,Arguments...)
%
%   SpecialParamIndex : index into the Arguments which denotes the (variable) regularization 
%                       pararameter. Can also be a vector of indices, if there are multiple such 
%                       parameters.
%
%   Arguments... : optional list of arguments to the LearningFunction (passed after the model)
%
% Out:
%   Model : the newly learned model
%
% Notes:
%   The global variable tracking.cache.max_cached_models can be used to control how many models may
%   be held in cache simultaneously (default: 20)
%
% Examples:
%   % suppose a learning function defined as follows exist:
%      newmodel = function mylearner(initialmodel,X,y,alpha,beta,gamma)
%
%   % ... which accepts an initial model as parameter, followed by several other parameters. 
%   % The Parameter "alpha" shall be special, and denote the regularization parameter (i.e., related 
%   % problems only differ in the value assigned to this parameter, whereas changes to any other 
%   % parameter denote an unrelated problem).
%
%   % Then, instead of invoking the function as:
%    result1 = mylearner([],myX,myY,0.1,mybeta,mygamma);
%    result2 = mylearner(result1,myX,myY,0.2,mybeta,mygamma);
%    result3 = mylearner(result2,myX,myY,0.3,mybeta,mygamma);
%   % It can be invoked as:
%    result1 = utl_warmstart([],@mylearner,3,myX,myY,0.1,mybeta,mygamma);
%    result2 = utl_warmstart([],@mylearner,3,myX,myY,0.2,mybeta,mygamma);
%    result3 = utl_warmstart([],@mylearner,3,myX,myY,0.3,mybeta,mygamma);
% 
%   % ... for any setting of its parameters - however, whenever the values of myX, myY, mybeta and 
%   % mygamma match the values of some previous invocation (perhaps with a different myalpha), the 
%   % model that was produced during that computation is being used instead of initial_model (and
%   % saved for the next case).
% 
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-10


persistent previous_models; % struct holding previous models (indexed by problem id / hash)
persistent previous_times;  % struct holding model creation times (indexed like previous_models)

% find a hash/id for the given learning problem
indices = 1:length(varargin); indices(parameter_index) = [];
hash = hlp_fingerprint([{learning_function},varargin(indices)]);
problem_id = sprintf('x%.0f',hash);

% see if we have a previous model stored for it
if isfield(previous_models,problem_id)    
    try
        % yes: try to use it
        model = learning_function(previous_models.(problem_id),varargin{:});
    catch e
        % no: fall back to the initial model in case of an error
        disp('Error evaluation learning function with cached model; falling back to initial model. Traceback:');
        env_handleerror(e);
        model = learning_function(initial_model,varargin{:});
    end
else
    % use the initial model right away
    model = learning_function(initial_model,varargin{:});
end

% store the model for the next time
previous_models.(problem_id) = model;
previous_times.(problem_id) = cputime;

global tracking;
if ~isfield(tracking.cache,'max_cached_models')
    tracking.cache.max_cached_models = 20; end

problems = fieldnames(previous_times);
% too many models in cache?
if length(problems) > tracking.cache.max_cached_models
    % get corresponding creation times
    times = struct2cell(previous_times);
    % sort them from oldest to newest
    [sorted,indices] = sort([times{:}]); %#ok<ASGLU>
    % get the indices of those that need to be removed
    toremove = indices(1:(length(problems) - tracking.cache.max_cached_models));
    % remove them
    previous_models = rmfield(previous_models,problems(toremove));
    previous_times = rmfield(previous_times,problems(toremove));
end
