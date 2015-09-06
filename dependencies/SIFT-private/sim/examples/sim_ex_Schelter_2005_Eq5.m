function [expr morder] = sim_ex_Schelter_2005_Eq5(varargin)
% Simulation:  Schelter 2005 Eq 5
%
% Description:  
% 
% 5-variate VAR[4] system of coupled oscillators. 
% This system was first described in [1].
% 
% The directed graph for this model is:
%
% x2 -> x1
% x3 -> x2
% x3 -> x4
% x3 -> x5
% x4 -> x2
% x5 -> x3
% x5 -> x4
%
% The dependency structure of this model can 
% be viewed by executing the following command:
%
% >>hlp_viewGraphicsResource('sim/Schelter_2005_Eq5.jpg');
%
% Author Credits:
% 
% Tim Mullen, 2011
%
% References and Code:
%
% [1] (Eq. 5) Schelter B, Winterhalder M, Eichler M, Peifer M, Hellwig B, Guschlbauer B, Lucking CH, Dahlhaus R, Timmer J (2005) Testing for directed influences among neural signals using partial directed coherence. Journal of neuroscience methods 152:210-9
%
% ------------------------------------------------------------------------

% specify the default system of equations
expr_def = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+  e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2) + -0.3*x3(t-4) + 0.6*x4(t-1) + e2(t)' ...
    'x3(t) = 0.8*x3(t-1) + -0.7*x3(t-2) + -0.1*x5(t-3) + e3(t)' ...
    'x4(t) = 0.5*x4(t-1) +  0.9*x3(t-2) +  0.4*x5(t-2) + e4(t)' ...
    'x5(t) = 0.7*x5(t-1) + -0.5*x5(t-2) + -0.2*x3(t-1) + e5(t)' ...
};

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},4,[1 Inf],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end