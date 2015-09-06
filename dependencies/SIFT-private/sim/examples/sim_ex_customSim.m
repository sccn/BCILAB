function [expr morder] = sim_ex_customSim(varargin)
% Simulation:  Custom Simulation
%
% Description:  
% 
% Custom simulation
%
% Recommended Settings:
%
% Author Credits:
% 
% N/A
%
% References and Code:
%
% N/A
%
% ------------------------------------------------------------------------

% specify the default system of equations
expr_def = {''};

% set up argument definitions
arg_define(varargin, ...
    arg({'expr','DynamicalEquations'},expr_def,[],'System of equations','type','cellstr','shape','row'), ...
    arg_subtoggle({'savesim','SaveSimulation'},'off', ...
    { ...
        arg({'simName','SimName'},'myCustomSim','','Simulation name. This is a user-friendly name for the simulation (selectable in GUI). This will also be used to auto-generate an m-file for the simulation. Simulation will be saved as /sim/examples/sim_ex_[SimNameSafe].m  where SimNameSafe is a filename-safe version of SimName'), ...
        arg({'simDescrip','SimDescription'},'','','Description for the simulation','type','char','shape','row'), ...
        arg({'simAuthor','SimAuthor'},'N/A','','Author credits for simulation'), ...
        arg({'simRef','SimRef'},'N/A','','Reference(s) for the simulation'), ... 
    },'Save the simulation'), ...
    arg({'morder','ModelOrder'},[],[],'Model order. This is mandatory'));

if isempty(morder)
    error('SIFT:sim_examples:badParam','ModelOrder must be specified. Make sure this matches the order of the dynamical equations');
end

    