function EEG = hlp_laplaceRef(EEG,varargin)
% Perform laplacian re-referecing using the spherical spline method of Perrin et al 1989
% This computes scalp current density or current source density
%
% This function requires the CSD toolbox by Jorgen Kayser
% (http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox)
%
% Author: Tim Mullen, SCCN/INC/UCSD Dec, 2013

% check if CSD toolbox exists
if ~exist('CSD','file') || ~exist('GetGH.m','file')
    error('CSD toolbox not installed.  Please download from [http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox] and place in Matlab path');
end

% parse inputs
g = finputcheck(varargin, ...
            {'sphSplineFlex' 'real' {} 4    ...
             'headRadius',   'real' {} 10   ...
             'smoothLambda', 'real' {} 1.0e-5 ...
            },'ignore','quiet');

% set up paths
locsdir = [tempdir() 'tmplocs.locs'];
csddir  = [tempdir() 'tmplocs.csd'];

% write the chanlocs to a temp location
writelocs(EEG.chanlocs,locsdir,'filetype','loc');

% convert chanlocs to CSD toolbox format
evalstr = sprintf(['ConvertLocations(' hlp_tostring(locsdir) ',' hlp_tostring(csddir) ');']);
evalc(evalstr);

% read locations
[~,M] = evalc(['ExtractMontage(' hlp_tostring(csddir) ',' hlp_tostring({EEG.chanlocs.labels}') ');']);

% clean up
delete(locsdir);
delete(csddir);

% make G and H matrices
[G H] = GetGH(M,g.sphSplineFlex);

% compute current source density
EEG.data = double(CSD(single(EEG.data),G,H));

