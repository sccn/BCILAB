%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to simulate various time-varying MVAR processes              %%%
%%% Author: Tim Mullen, 2011, SCCN, INC, UCSD                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear the workspace
clear;

%% STEP 1: Define the parameters for the system of coupled oscillators

%% Example 1: Simple bivariate coupled oscillator example

SamplingRate = 100;      % Hz
Nl = 50;                % length of each epoch (samples)
Nr = 50;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 2;          % VAR model order
 
expr = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+ e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2) + e2(t)' ...
};

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);

%% Example 1b: trivariate coupled oscillators with non-stationary (1-Hz sinusoidal) coupling dynamics


SamplingRate = 100;      % Hz
Nl = 500;                % length of each epoch (samples)
Nr = 100;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 2;          % VAR model order
% ModelOrder = 4;        % Model Order we will use
f0 = 10;                 % central oscillation frequency
expr = {...
    ['x1(t) = ' sim_dampedOscillator(f0,10,100,1) '                                    + e1(t)'] ... 
    ['x2(t) = ' sim_dampedOscillator(f0,2,100,2) ' + -0.1*x1(t-2)                      + e2(t)'] ...
    ['x3(t) = ' sim_dampedOscillator(f0,2,100,3) ' + {0.3*sin(2*pi*t/100)+0.3}*x1(t-2) + e3(t)'] ...
};

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);


%% Example 2:  System of coupled oscillators 
% from (Ex 3.1, eq 11-15) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed 
% influences among neural signals using renormalized partial directed coherence. 
% Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

SamplingRate = 100;      % Hz

Nl = 500;                % length of each epoch (samples)
Nr = 100;                % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 3;          % VAR model order

% write the system of equations for coupled damped oscillators
expr = { ...
    'x1(t) = 0.9*x1(t-1)  + 0.3*x2(t-2)  + e1(t)' ...
    'x2(t) = 1.3*x2(t-1)  + -0.8*x2(t-2) + e2(t)' ...
    'x3(t) = 0.3*x1(t-2)  + 0.6*x2(t-1)  + e3(t)' ...
    'x4(t) = -0.7*x4(t-3) + -0.7*x1(t-3) + 0.3*x5(t-3) + e4(t)' ...
    'x5(t) = 1*x5(t-1)    + -0.4*x5(t-2) + 0.3*x4(t-2) + e5(t)' ...
    };
  
% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);


%% Example 3: Static system of coupled oscillators 
% (Ex 3.2, eq 16-20) Schelter B, Timmer J, Eichler M (2009) Assessing the strength of directed 
% influences among neural signals using renormalized partial directed coherence. 
% Journal of neuroscience methods 179:121-30 Available at: http://www.ncbi.nlm.nih.gov/pubmed/19428518.

SamplingRate = 10;          % Hz

Nl = 500;                   % length of each epoch (samples)
Nr = 100;                   % number of trials (realisations)
ndisc = 1000;               % number of samples to discard from VAR model (startup transients)
ModelOrder = 2;         % VAR model order

% write the system of equations for coupled damped oscillators
expr = { ...
    'x1(t) = 1.9*x1(t-1) + -0.999*x1(t-2) + e1(t)' ...
    'x2(t) = 0.9*x2(t-2) + -0.2*x1(t-1)   + e2(t)' ...
    'x3(t) = -0.3*x3(t-1)+ 0.4*x4(t-1)    + -0.3*x5(t-2) + e3(t)' ...
    'x4(t) = 1.3*x4(t-1) + -0.7*x4(t-2)   + e4(t)' ...
    'x5(t) = 0.7*x5(t-2) + 0.3*x1(t-1)    + e5(t)' ...
    };

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);


%% Example 4: Static system of coupled oscillators
% (eq 5) from Schelter B, Winterhalder M, Eichler M, Peifer M, Hellwig B, Guschlbauer B, 
% Lucking CH, Dahlhaus R, Timmer J (2005) Testing for directed influences among neural signals using 
% partial directed coherence. Journal of neuroscience methods 152:210-9 
% Available at: http://www.ncbi.nlm.nih.gov/pubmed/16269188.

SamplingRate = 100;         % Hz
Nl = 100;                   % length of each epoch (samples)
Nr = 10;                   % number of trials (realisations)
ndisc = 1000;               % number of samples to discard from VAR model (startup transients)
ModelOrder = 4;             % VAR model order

% write the system of equations for coupled damped oscillators
expr = {...
    'x1(t) = 0.6*x1(t-1) +  0.65*x2(t-2)+e1(t)' ... 
    'x2(t) = 0.5*x2(t-1) + -0.3*x2(t-2)+ -0.3*x3(t-4) +0.6*x4(t-1) + e2(t)' ...
    'x3(t) = 0.8*x3(t-1) + -0.7*x3(t-2)+ -0.1*x5(t-3) + e3(t)' ...
    'x4(t) = 0.5*x4(t-1) +  0.9*x3(t-2)+ 0.4*x5(t-2)  + e4(t)' ...
    'x5(t) = 0.7*x5(t-1) + -0.5*x5(t-2)+ -0.2*x3(t-1) + e5(t)' ...
};

% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);

%% Example 5: Single-trial, Time-Varying Simulation (Mullen, 2011 unpublished)
%  System of stocastically-forced, damped coupled oscillators with
%  time-varying coupling coefficients
%  This simulation creates a simualted "seizure" with time-varying coupling
%  between clusters of sources which switches between 4 different stages.

SamplingRate = 100; % Hz

Nl = 5*60*SamplingRate;  % length of each epoch (5 minutes)
Nr = 1;                  % number of trials (realisations)
ndisc = 1000;            % number of samples to discard from VAR model (startup transients)
ModelOrder = 6;          % VAR model order

% Set the fundamental frequency (f) and damping time (tau) of the 
% oscillators for each cluster of sources
% NOTE: A small tau generates a "noisier" signal

% Cluster 1
f1=20; 
tau1 = 20;

% Cluster 2
f2=20;
tau2 = 3;

% Cluster 3
f3=10;
tau3 = 7;

% Cluster 4
f4 = 10;
tau4 = 6;


% set the approximate durations of each stage
S1_width = 5*SamplingRate;
S2_width = 5*SamplingRate;
S3_width = 5*SamplingRate;
S4_width = 5*SamplingRate;

S0 = 60*SamplingRate; % start of seizure (1 minute)

S1_center = S0+S1_width/2;                      % center of stage 1
S2_center = S1_center+S1_width;                 % center of stage 2
S3_center = S2_center;                          % center of stage 3
S4_center = S3_center+S3_width;                 % center of stage 4

Offset = 0;

% write the system of equations for coupled damped oscillators
expr = { ...
      sprintf('x1(t) = {2*exp(-1/(%f+normpdfg(t,%f,%f,%f,100)))*cos(2*pi*20.000000/100.000000)}*x1(t-1) + {-exp(-2/(%f+normpdfg(t,%f,%f,%f,100)))}*x1(t-2) + e1(t)',tau1,S1_width/2,8,Offset+S1_center, tau1,S1_width,8,Offset+S1_center) ...                    % Ictal driver
      sprintf('x2(t) = %s + {normpdfg(t,%f,%f,%f,0.01)}*x3(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,1.3)}*x1(t-6)    + e2(t)',sim_dampedOscillator(f2,tau2,SamplingRate,2),  S1_width/2,8,Offset+S1_center,   S1_width/2,8,Offset+S1_center,  S1_width/2,8,Offset+S1_center) ...                     % CLUSTER 2
      sprintf('x3(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x4(t-2)  + {normpdfg(t,%f,%f,%f,0.3)}*x5(t-3)    + e3(t)',sim_dampedOscillator(f2-1,tau2,SamplingRate,3),  S1_width/2,8,Offset+S1_center,   S1_width/2,8,Offset+S1_center,   (S2_width+S3_width)/2,10,Offset+S3_center) ...       % CLUSTER 2
      sprintf('x4(t) = %s +  {normpdfg(t,%f,%f,%f,0.7)}*x2(t-2) + {normpdfg(t,%f,%f,%f,0.1)}*x3(t-2)                                          + e4(t)',sim_dampedOscillator(f2+1,tau2,SamplingRate,4),  S1_width/2,8,Offset+S1_center,  S1_width/2,8,Offset+S1_center) ... % CLUSTER 2
      sprintf('x5(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x6(t-2) + {normpdfg(t,%f,%f,%f,0.5)}*x3(t-3)        + e5(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,5), S3_width/2,8,Offset+S3_center , S3_width/2,8,Offset+S3_center) ...  % CLUSTER 3
      sprintf('x6(t) = %s +  {normpdfg(t,%f,%f,%f,0.001)}*x5(t-2)                                             + e6(t)'  ,sim_dampedOscillator(f3,tau3,SamplingRate,6), S3_width/2,8,Offset+S3_center) ...                                  % CLUSTER 3
      sprintf('x7(t) = %s + {normpdfg(t,%f,%f,%f,1.3)}*x4(t-6) + {normpdfg(t,%f,%f,%f,0.01)}*x9(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x8(t-2) + {normpdfg(t,%f,%f,%f,0.01)}*x10(t-2)    + e7(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,7),     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center,     S4_width/2,8,Offset+S4_center) ...  % CLUSTER 4
      sprintf('x8(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e8(t)' ,sim_dampedOscillator(f4,tau4,SamplingRate,8),     S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x9(t) = %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e9(t)' ,sim_dampedOscillator(f4+1,tau4,SamplingRate,9),   S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x10(t)= %s + {normpdfg(t,%f,%f,%f,0.6)}*x7(t-2)         + e10(t)',sim_dampedOscillator(f4-1,tau4,SamplingRate,10),  S4_width/2,8,Offset+S4_center) ... % CLUSTER 4
      sprintf('x11(t)= 1.3*x11(t-1) + -0.8*x11(t-2)                             + e11(t)') ...    % Decoupled
      sprintf('x12(t)= 1.2*x12(t-1) + -0.4*x12(t-2)                             + e12(t)') ...    % Decoupled
      sprintf('x13(t)= 0.8*x13(t-1) + -0.4*x13(t-2) + -0.1*x13(t-4)             + e13(t)') ...    % Decoupled
      };
  
% create prototype VAR structure
Aproto = sim_genVARModelFromEq(expr,ModelOrder);





%% STEP 2: Simulate the VAR process

[A] = sim_genTVARcoeffs(Aproto,ModelOrder,Nl,'NumSamplesToDiscard',ndisc,'Verbose',true);

%% STEP 3: generate data from the VAR model

% Specify the noise covariance matrix. 
% Sigma is the noise variance.
sigma = 1;
M = size(A{1},1);
C = sigma*eye(M);             

% % hyperbolic secant noise
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2/pi,0,'hsec'),[2 1 3]);

% laplacian noise (generalized gaussian)
% data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,1,1,'gengauss'),[2 1 3]);

% gaussian noise
data = permute(tvarsim(zeros(1,M),A,C,[Nl Nr],ndisc,2,1,'gengauss'),[2 1 3]);


%% STEP 4: Create EEGLAB dataset and source potentials

EEG = eeg_emptyset;
EEG.data = data;
[EEG.icaweights EEG.icasphere EEG.icawinv] = deal(eye(M));
EEG.icaact = [];
EEG.srate = SamplingRate;
EEG.times = ((0:(Nl-1))/SamplingRate)*1000;   % ms
EEG.pnts = Nl;
EEG.trials = Nr;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end)/1000;  % sec
EEG.nbchan = M;
EEG.setname = 'VAR Simulation';
EEG.condition = 'VAR Simulation';

%%
EEG = eeg_checkset(EEG);
pop_eegplot(EEG);

ALLEEG = EEG;
CURRENTSET = length(ALLEEG);
eeglab redraw;

%% STEP 6: Proceed through analysis and visualization as directed by Chapter 6 of the manual 
% (or see ScriptingExample.m for command-line alternatives)