function res = exp_lambda(body,varargin)
% Create an expression from an anonymous function (with optional arguments).
% Lambda-Expression = exp_lambda(Function-Body, Arguments)

res = struct('head',{body},'parts',{varargin});