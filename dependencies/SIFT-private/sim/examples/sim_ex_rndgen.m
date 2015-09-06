function [expr morder] = sim_ex_rndgen(varargin)
% Simulation:  Random Model Generation
%
% Description:  
% 
% Generate a stable, random VAR model
%
% Recommended Settings:
% 
% For large models, it is recommended that you keep the model quite sparse
% (< 1% interacting variables) otherwise searching for a stable model may
% take a long time.
%
% Author Credits:
% 
% Tim Mullen
%
% References and Code:
%
% See sim_genRndVARcoeffs()
%
% ------------------------------------------------------------------------

expr    = {};
morder  = [];
% set up argument definitions
arg_define(varargin,arg_sub({'rndopts','Options'},{},@sim_genRndVARcoeffs,'Options for random model generation'));


    
