function [expr morder] = sim_ex_Schelter_2009_3_1(varargin)
% Simulation:  Schelter 2009 Eq 3.1
%
% Description:  
% 
% 5-variate VAR[3] system of coupled oscillators. 
% This system was first described in [1].
% 
% The directed graph for this model is:
% x1 -> x3
% x1 -> x4
% x2 -> x1
% x2 -> x3
% x4 -> x5
% x5 -> x4
%
% The dependency structure of this model can 
% be viewed by executing the following command:
%
% >>hlp_viewGraphicsResource('sim/Schelter_2009_3_1.jpg');
%
% Author Credits:
% 
% Tim Mullen, 2011
%
% References and Code:
%
% [1] (Ex 3.1, Eq. 11-15) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed influences among neural signals using renormalized partial directed coherence. Journal of neuroscience methods 179:121-30
%
% ------------------------------------------------------------------------

% specify the default system of equations
expr_def = { ...
    'x1(t) = 0.9*x1(t-1)  + 0.3*x2(t-2)  + e1(t)' ...
    'x2(t) = 1.3*x2(t-1)  + -0.8*x2(t-2) + e2(t)' ...
    'x3(t) = 0.3*x1(t-2)  + 0.6*x2(t-1)  + e3(t)' ...
    'x4(t) = -0.7*x4(t-3) + -0.7*x1(t-3) + 0.3*x5(t-3) + e4(t)' ...
    'x5(t) = 1*x5(t-1)    + -0.4*x5(t-2) + 0.3*x4(t-2) + e5(t)' ...
    };

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},3,[1 Inf],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end