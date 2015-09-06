function [expr morder] = sim_ex_Schelter_2009_3_2(varargin)
% Simulation:  Schelter 2009 Eq 3.2
%
% Description:  
% 
% 5-variate VAR[2] system of coupled oscillators. 
% This system was first described in [1].
% 
% The directed graph for this model is:
% x1 -> x2
% x1 -> x5
% x4 -> x3
% x5 -> x3
%
% The dependency structure of this model can 
% be viewed by executing the following command:
%
% >>hlp_viewGraphicsResource('sim/Schelter_2009_3_2.jpg');
%
% Author Credits:
% 
% Tim Mullen, 2011
%
% References and Code:
%
% [1] (Ex 3.2, Eq. 16-20) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed influences among neural signals using renormalized partial directed coherence. Journal of neuroscience methods 179:121-30
%
% ------------------------------------------------------------------------

% specify the default system of equations
expr_def = { ...
    'x1(t) =  1.9*x1(t-1) + -0.999*x1(t-2) + e1(t)' ...
    'x2(t) =  0.9*x2(t-2) + -0.2*x1(t-1)   + e2(t)' ...
    'x3(t) = -0.3*x3(t-1) +  0.4*x4(t-1)   + -0.3*x5(t-2) + e3(t)' ...
    'x4(t) =  1.3*x4(t-1) + -0.7*x4(t-2)   + e4(t)' ...
    'x5(t) =  0.7*x5(t-2) +  0.3*x1(t-1)   + e5(t)' ...
    };

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},2,[1 Inf],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end