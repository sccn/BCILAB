function [expr morder] = sim_ex_triCoupledOscillator(varargin)
% Simulation:  Trivariate Coupled Oscillator
%
% Description:  
% 
% Trivariate 10-Hz coupled oscillators with non-
% stationary (1-Hz sinusoidal) coupling dynamics
%
% Recommended Settings:
% Sampling Rate: 100 Hz
%
% The directed graph for this model can be viewed by executing the following command:
%
% >>hlp_viewGraphicsResource('sim/triCoupledOscillator.jpg');
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
f0 = 10;
expr_def = {...
    ['x1(t) = ' sim_dampedOscillator(f0,10,100,1) '                                    + e1(t)'] ... 
    ['x2(t) = ' sim_dampedOscillator(f0,2,100,2) ' + -0.1*x1(t-2)                      + e2(t)'] ...
    ['x3(t) = ' sim_dampedOscillator(f0,2,100,3) ' + {0.3*sin(2*pi*t/100)+0.3}*x1(t-2)   + e3(t)'] ...
};

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},2,[],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end
