% startupadigator: will add adigator home directory and util directory to
% the matlab path
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
currpath = pwd;
addpath([currpath,filesep,'lib']);
addpath([currpath,filesep,'util']);

fprintf('ADiGator Successfully Installed\n')