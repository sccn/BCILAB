function init_unlocbox()
%INIT_TOOLBOX Initialize the toolbox
%   Usage: init_unlocbox()
%
%   Initialisation script for the convex optimization problems toolbox
%   This script add the different path needed to run the toolbox
%
%   Url: http://unlocbox.sourceforge.net/doc/init_unlocbox.php

% Copyright (C) 2012 LTS2-EPFL, by Nathanael Perraudin, Gilles Puy,
% David Shuman, Pierre Vandergheynst.
% This file is part of UnLocBoX version 1.1.70
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Author: Nathanael Perraudin
% E-mail: nathanael.perraudin@epfl.ch
% Date: nov 2012


%% adding dependency
global GLOBAL_useGPU;
global GLOBAL_path;
GLOBAL_path = fileparts(mfilename('fullpath'));
GLOBAL_useGPU = 0;

addpath(genpath(GLOBAL_path));

% Load the version number
bp=[GLOBAL_path,filesep];
[FID, MSG] = fopen ([bp,'unlocbox_version'],'r');
if FID == -1
    error(MSG);
else
    unlocbox_version = fgetl (FID);
    fclose(FID);
end

banner = sprintf(strcat(... 
'UnLocBoX version %s. Copyright 2012-2013 LTS2-EPFL, by Nathanael Perraudin'), ...
                   unlocbox_version);
  

if GLOBAL_useGPU
    init_gpu();
end

% display banner
disp(banner);

