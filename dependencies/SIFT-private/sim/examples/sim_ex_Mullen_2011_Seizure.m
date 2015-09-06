function [expr morder] = sim_ex_Mullen_2011_Seizure(varargin)
% Simulation:  Epileptic Seizure
%
% Description:  
% 
%   13-variate VAR[6] system of stocastically-forced, damped coupled oscillators with time-varying (non-stationary) coupling.
%
%   This simulation creates a simulated "seizure" with time-varying coupling between clusters of sources which switches between 4 different stages.
%
%   The simulation is designed to be single-trial.
% 
%   Recommended Settings: 
%   Number of trials: 1
%   Sampling Rate:    100
%
%   The directed graph for this model can be viewed by executing the following command:
%
%   >>hlp_viewGraphicsResource('/sim/Mullen_2011_Seizure.jpg');
%
%
% Author Credits:
% 
%   Tim Mullen, 2011
%
% References and Code:
%
%   Mullen, 2011. Unpublished
%
% ------------------------------------------------------------------------

% Cluster 1
f1   = .2;      % central frequency
tau1 = .2;      % damping time (larger-->highest variance)
% Cluster 2
f2   = .20;
tau2 = 3;
% Cluster 3
f3   = .10;
tau3 = 7;
% Cluster 4
f4   = .10;
tau4 = 6;

% start of seizure (1 minute)
S0 = 6000; 
% set the approximate durations of each stage
S1_duration = 500;  % 5 sec at 100 Hz
S2_duration = 500;
S3_duration = 500;
S4_duration = 500;
% set seizure stage mid-points
S1_center = S0+S1_duration/2;                      % center of stage 1
S2_center = S1_center+S1_duration;                 % center of stage 2
S3_center = S2_center;                             % center of stage 3
S4_center = S3_center+S3_duration;                 % center of stage 4

Offset = 0;
SamplingRate = 1;

% specify the default system of equations
expr_def = { ...
          sprintf('x1(t) = {2*exp(-1/(%f+normpdfg(t,%f,%f,%f,100)))*cos(2*pi*%f/1)}*x1(t-1) + {-exp(-2/(%f+normpdfg(t,%f,%f,%f,100)))}*x1(t-2) + e1(t)',tau1,S1_duration/2,8,Offset+S1_center, f1, tau1,S1_duration,8,Offset+S1_center) ...                    % Ictal driver
          sprintf('x2(t) = %s + {normpdfg(t,%f,%f,%f,0.01)}*x3(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,1.3)}*x1(t-6)    + e2(t)',sim_dampedOscillator(f2,tau2,SamplingRate,2),  S1_duration/2,8,Offset+S1_center,   S1_duration/2,8,Offset+S1_center,  S1_duration/2,8,Offset+S1_center) ...                     % CLUSTER 2
          sprintf('x3(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,0.3)}*x5(t-3)    + e3(t)',sim_dampedOscillator(f2-1,tau2,SamplingRate,3),  S1_duration/2,8,Offset+S1_center,   S1_duration/2,8,Offset+S1_center,   (S2_duration+S3_duration)/2,10,Offset+S3_center) ...       % CLUSTER 2
          sprintf('x4(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x3(t-2)                                          + e4(t)',sim_dampedOscillator(f2+1,tau2,SamplingRate,4),  S1_duration/2,8,Offset+S1_center,  S1_duration/2,8,Offset+S1_center) ... % CLUSTER 2
          sprintf('x5(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x6(t-2) + {normpdfg(t,%f,%f,%f,0.5)}*x3(t-3)        + e5(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,5), S3_duration/2,8,Offset+S3_center , S3_duration/2,8,Offset+S3_center) ...  % CLUSTER 3
          sprintf('x6(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x5(t-2)                                             + e6(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,6), S3_duration/2,8,Offset+S3_center) ...                                  % CLUSTER 3
          sprintf('x7(t) = %s + {normpdfg(t,%f,%f,%f,1.3)}*x4(t-6) + {normpdfg(t,%f,%f,%f,0.01)}*x9(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x8(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x10(t-2)    + e7(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,7),     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center) ...  % CLUSTER 4
          sprintf('x8(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e8(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,8),     S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x9(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e9(t)' ,sim_dampedOscillator(f4+1,tau4,SamplingRate,9),   S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x10(t)= %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e10(t)',sim_dampedOscillator(f4-1,tau4,SamplingRate,10),  S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x11(t)= 1.3*x11(t-1) + -0.8*x11(t-2)                             + e11(t)') ...    % Decoupled
          sprintf('x12(t)= 1.2*x12(t-1) + -0.4*x12(t-2)                             + e12(t)') ...    % Decoupled
          sprintf('x13(t)= 0.8*x13(t-1) + -0.4*x13(t-2) + -0.1*x13(t-4)             + e13(t)') ...    % Decoupled
          };

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations'), ...
    arg({'morder','ModelOrder'},6,[1 Inf],'Model order. This is mandatory'), ...
    arg_subtoggle({'setDynamics','SetDynamics'},'off', ...
        {...
            arg_sub({'c1','Cluster1'},{},...
            { ...
                arg({'f0','Freq'},f1,[0 Inf],'Central frequency (Hz). This is in normalized Hz (i.e. f0=0.5=SamplingRate/2).'), ...
                arg({'tau','DampingTime'},tau1,[0 Inf],'Damping time') ...
            },'Cluster 1 Settings'), ...
            arg_sub({'c2','Cluster2'},{},...
            { ...
                arg({'f0','Freq'},f2,[0 Inf],'Central frequency (Hz). This is in normalized Hz (i.e. f0=0.5=SamplingRate/2).'), ...
                arg({'tau','DampingTime'},tau2,[0 Inf],'Damping time') ...
            },'Cluster 2 Settings'), ...
            arg_sub({'c3','Cluster3'},{},...
            { ...
                arg({'f0','Freq'},f3,[0 Inf],'Central frequency (Hz). This is in normalized Hz (i.e. f0=0.5=SamplingRate/2).'), ...
                arg({'tau','DampingTime'},tau3,[0 Inf],'Damping time') ...
            },'Cluster 3 Settings'), ...
            arg_sub({'c4','Cluster4'},{},...
            { ...
                arg({'f0','Freq'},f4,[0 Inf],'Central frequency (Hz). This is in normalized Hz (i.e. f0=0.5=SamplingRate/2).'), ...
                arg({'tau','DampingTime'},tau4,[0 Inf],'Damping time') ...
            },'Cluster 4 Settings'), ...
            arg_sub({'s1','Stage1'},{},...
            { ...
                arg({'duration','Duration'},S1_duration,[0 Inf],'Duration (samples)'), ...
                arg({'center','Center'},S1_center,[0 Inf],'Midpoint of stage (samples)'), ...
            },'Stage 1 Settings'), ...
            arg_sub({'s2','Stage2'},{},...
            { ...
                arg({'duration','Duration'},S2_duration,[0 Inf],'Duration (samples)'), ...
                arg({'center','Center'},S2_center,[0 Inf],'Midpoint of stage (samples)'), ...
            },'Stage 2 Settings'), ...
            arg_sub({'s3','Stage3'},{},...
            { ...
                arg({'duration','Duration'},S3_duration,[0 Inf],'Duration (samples)'), ...
                arg({'center','Center'},S3_center,[0 Inf],'Midpoint of stage (samples)'), ...
            },'Stage 3 Settings'), ...
            arg_sub({'s4','Stage4'},{},...
            { ...
                arg({'duration','Duration'},S4_duration,[0 Inf],'Duration (samples)'), ...
                arg({'center','Center'},S4_center,[0 Inf],'Midpoint of stage (samples)'), ...
            },'Stage 4 Settings'), ...
            arg({'szStart','SeizureStart'},S0,[0 Inf],'Start of seizure (samples)'), ...
        },'Override default system dynamics'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified');
end

if g.setDynamics.arg_selection
    % re-evaluate expression with user-defined settings
    
    % set central frequency and damping times
    % Cluster 1
    f1   = g.setDynamics.c1.f0;      % central frequency
    tau1 = g.setDyamics.c1.tau;      % damping time (larger-->highest variance)
    % Cluster 2
    f2   = g.setDynamics.c2.f0;
    tau2 = g.setDyamics.c2.tau;
    % Cluster 3
    f3   = g.setDyamics.c3.f0;
    tau3 = g.setDyamics.c3.tau;
    % Cluster 4
    f4   = g.setDyamics.c4.f0;
    tau4 = g.setDyamics.c4.tau;
    
    % set start of seizure
    S0 = g.setDynamics.szStart;
    
    % set the approximate durations of each stage
    S1_duration = g.setDynamics.s1.duration;
    S2_duration = g.setDynamics.s2.duration;
    S3_duration = g.setDynamics.s3.duration;
    S4_duration = g.setDynamics.s4.duration;
    
    % set seizure stage mid-points
    if ~isempty(g.setDynamics.s1.center)
        S1_center = g.setDynamics.s1.center;
    else
        S1_center = S0+S1_duration/2;      
    end
    if ~isempty(g.setDynamics.s2.center)
        S2_center = g.setDynamics.s2.center;
    else
        S2_center = S1_center+S1_duration;    
    end
    if ~isempty(g.setDynamics.s3.center)
        S3_center = g.setDynamics.s3.center;
    else
         S3_center = S2_center;      
    end
    if ~isempty(g.setDynamics.s4.center)
        S4_center = g.setDynamics.s4.center;
    else
        S4_center = S3_center+S3_duration;     
    end
       
    
    % generate system of equations
    expr = { ...
          sprintf('x1(t) = {2*exp(-1/(%f+normpdfg(t,%f,%f,%f,100)))*cos(2*pi*%f/1)}*x1(t-1) + {-exp(-2/(%f+normpdfg(t,%f,%f,%f,100)))}*x1(t-2) + e1(t)',tau1,S1_duration/2,8,Offset+S1_center, f1, tau1,S1_duration,8,Offset+S1_center) ...                    % Ictal driver
          sprintf('x2(t) = %s + {normpdfg(t,%f,%f,%f,0.01)}*x3(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,1.3)}*x1(t-6)    + e2(t)',sim_dampedOscillator(f2,tau2,SamplingRate,2),  S1_duration/2,8,Offset+S1_center,   S1_duration/2,8,Offset+S1_center,  S1_duration/2,8,Offset+S1_center) ...                     % CLUSTER 2
          sprintf('x3(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,0.3)}*x5(t-3)    + e3(t)',sim_dampedOscillator(f2-1,tau2,SamplingRate,3),  S1_duration/2,8,Offset+S1_center,   S1_duration/2,8,Offset+S1_center,   (S2_duration+S3_duration)/2,10,Offset+S3_center) ...       % CLUSTER 2
          sprintf('x4(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x3(t-2)                                          + e4(t)',sim_dampedOscillator(f2+1,tau2,SamplingRate,4),  S1_duration/2,8,Offset+S1_center,  S1_duration/2,8,Offset+S1_center) ... % CLUSTER 2
          sprintf('x5(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x6(t-2) + {normpdfg(t,%f,%f,%f,0.5)}*x3(t-3)        + e5(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,5), S3_duration/2,8,Offset+S3_center , S3_duration/2,8,Offset+S3_center) ...  % CLUSTER 3
          sprintf('x6(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x5(t-2)                                             + e6(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,6), S3_duration/2,8,Offset+S3_center) ...                                  % CLUSTER 3
          sprintf('x7(t) = %s + {normpdfg(t,%f,%f,%f,1.3)}*x4(t-6) + {normpdfg(t,%f,%f,%f,0.01)}*x9(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x8(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x10(t-2)    + e7(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,7),     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center,     S4_duration/2,8,Offset+S4_center) ...  % CLUSTER 4
          sprintf('x8(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e8(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,8),     S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x9(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e9(t)' ,sim_dampedOscillator(f4+1,tau4,SamplingRate,9),   S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x10(t)= %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e10(t)',sim_dampedOscillator(f4-1,tau4,SamplingRate,10),  S4_duration/2,8,Offset+S4_center) ... % CLUSTER 4
          sprintf('x11(t)= 1.3*x11(t-1) + -0.8*x11(t-2)                             + e11(t)') ...    % Decoupled
          sprintf('x12(t)= 1.2*x12(t-1) + -0.4*x12(t-2)                             + e12(t)') ...    % Decoupled
          sprintf('x13(t)= 0.8*x13(t-1) + -0.4*x13(t-2) + -0.1*x13(t-4)             + e13(t)') ...    % Decoupled
          };
end