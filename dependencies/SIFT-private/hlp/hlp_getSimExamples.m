function sims = hlp_getSimExamples(varargin)
% return a cell array of simulation names and associated function
% handles. The format is
% sims = {{name1 func_handle1} {name2 func_handle2} ... }
% if defaultNameOnly = true, then we only return name1 as a string
%
% hlp_getSimExamples('defNameOnly') returns only the name of the 
%  default algorithm (as a string)
%
% hlp_getSimExamples('mfileNameOnly',simName) where simName is a string
%   with the human-readable name of a valid algorithm, returns the m-file 
%   name of the function implementing the specified approach
% 
% NOTES:
%
% If the preamble text of the function (c.f. hlp_getFcnPreambleText())
% begins with the line
%
% Simulation: <simulation_name>
%
% then <simulation_name> is returned as the human-readable simulation name. 
% Otherwise, if this line cannot be found, the function name is returned.
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


persistent simEx;

defaultNameOnly = false;
simName = '';

if nargin==1 && strcmp(varargin{1},'defaultNameOnly')
    defaultNameOnly = true;
end

if nargin==2 && ismember_bc('mfileNameOnly',varargin);
    simName = varargin{2};
    if ~ischar(simName)
        error('SIFT:sim_getSimExamples:badInput', ...
              'Bad argument pair for ''mfileNameOnly'''); 
    end
end

if ~isempty(simEx)
    % only do the check once (it's time-consuming)
    % ... retrive "cached" version
    sims = simEx;
else
    
    sims = {};
        
    % get the names of all mvar_* functions in the /est/mvar folder
    siftroot = fileparts(which('StartSIFT'));
    fpath    = [siftroot filesep 'sim' filesep 'examples' filesep];
    simFcns = wildcardsearch(fpath,'*sim_ex_*.m',true,true);
    simFcns = regexprep(simFcns,['.*examples' filesep],'');
    simFcns = regexprep(simFcns,'\.m','');
    
    % cycle through the list of algorithms
    for k=1:length(simFcns)
        
        % get the help text (H1) for the algorithm entry function
        try preText = hlp_getFcnPreambleText(simFcns{k});
        catch err
            disp(err.message);
            continue;
        end
        
        % extract the human-readable algorithm name from the text 
        % following the 'Simulation:' header
        preText = strtrim(preText);
        simName = regexpi(preText,'Simulation\s*[:]?\s*([^\n]*)','tokens');
        if ~isempty(simName)
            simName = strtrim(simName{1}{1});
        else
            simName = simFcns{k};
        end
        
        sims{end+1} = {simName str2func(simFcns{k})};
    end
    
    simEx = sims;
end

% return only the name of the first available algorithm
if defaultNameOnly
    %sims = 'Trivariate Coupled Oscillators';
    sims = sims{1}{1};
elseif ~isempty(simName)
    % get the function handle matching the desired algorithm name
    tmp = cellfun(@(x) fastif(isequal(x{1},simName),x{2},''),sims,'UniformOutput',false);
    tmp(cellfun(@isempty,tmp))=[];
    if isempty(tmp)
        error('SIFT:sim_getSimExamples:badSimulationName', ...
              'Unknown simulation ''%s''',simName);
    else
        % return function name as a string
        sims = char(tmp{1});
    end
end