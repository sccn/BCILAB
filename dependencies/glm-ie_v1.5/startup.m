% startup script to make Octave/Matlab aware of the glm-ie package
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 25

disp(['executing glm-ie startup script...']);

OCT = exist('OCTAVE_VERSION') ~= 0;           % check if we run Matlab or Octave

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located
if OCT && numel(mydir)==2 
  if strcmp(mydir,'./'), mydir = [pwd,mydir(2:end)]; end
end                 % OCTAVE 3.0.x relative, MATLAB and newer have absolute path

addpath(mydir(1:end-1))
addpath([mydir,'inf'])
addpath([mydir,'mat'])
if OCT, addpath([mydir,'mat/@matWav/private']), end
addpath([mydir,'mat'])
addpath([mydir,'pen'])
addpath([mydir,'pls'])
addpath([mydir,'pot'])
addpath([mydir,'doc'])

clear me mydir
