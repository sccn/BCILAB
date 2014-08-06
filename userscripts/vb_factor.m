function F = vb_factor(func,varargin)
% Define a factor in a Bayesian log posterior distribution.
% FactorDef = vb_factor(DistFunction,Arguments...)
%
% This function allows to define a factor that follows a given distribution and depends on parameters 
% of the posterior, possibly in a non-trivial way.
%
% In:
%   DistFunction : function handle for the respective distribution, e.g., @logmvnpdf
%
%   Arguments... : list of arguments to the distribution function, where arguments that depend
%                  on posterior parameters that shall be inferred are given as an anonymous 
%                  function that has a single argument with the name of the respective parameter.
%
% Out:
%   FactorDef : definition of the factor for use in a VB solver (e.g., vb_mvn)

F = struct('func',{func},'args',{varargin});
