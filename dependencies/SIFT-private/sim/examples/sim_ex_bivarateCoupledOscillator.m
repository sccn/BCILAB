function [expr morder] = sim_ex_bivarateCoupledOscillator(varargin)
% Simulation:  Bivariate Coupled Oscillator
%
% Description:  
% 
% Bivariate coupled oscillator example
%
% Author Credits:
% 
% Tim Mullen, 2011
%
% References and Code:
%
% N/A
%
% ------------------------------------------------------------------------

% specify the default system of equations
expr_def = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+ e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2) + e2(t)' ...
};

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},2,[1 Inf],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end
