function res = exp_hold(varargin)
% Hold the given expression from being evaluated; can be released with exp_releasehold.
% Held-Expression = Hold(Expression...)
%
% In:
%   Expression... : An expression or list of expressions
%
% Out:
%   Held-Expression ; The given expression(s), wrapped in a Hold expression
%
% See also:
%   exp_releasehold
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-19

res = struct('head',{@Hold},'parts',{varargin});