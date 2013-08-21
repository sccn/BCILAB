%
%
%  Setup_Nesta
% 
%  sets up the path for Nesta
%
%
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: May 2009
%
% NESTA NESTA Version 1.1
%   See also NESTA

% This adds a relative path, so it only works
% when you're in the NESTA directory

%addpath RecPF_v1.1/solver/
%addpath Misc/


% This adds an absolute path, so it works
% when the user is in any directory

p = mfilename('fullpath'); % the location of this file
path = fileparts(p);
addpath( fullfile(path,'RecPF_v1.1','solver' ) );
addpath( fullfile(path,'Misc') );
