function algs = hlp_getMVARalgorithms(varargin)
% return a cell array of MVAR algorithm names and associated function
% handles. The format is
% algs = {{name1 func_handle1} {name2 func_handle2} ... }
% if defaultNameOnly = true, then we only return name1 as a string
%
% hlp_getMVARalgorithms('defNameOnly') returns only the name of the 
%  default algorithm (as a string)
%
% hlp_getMVARalgorithms('mfileNameOnly',algName) where algName is a string
%   with the human-readable name of a valid algorithm, returns the m-file 
%   name of the function implementing the specified approach
% 
% NOTES:
%
% If the preamble text of the function (c.f. hlp_getFcnPreambleText())
% begins with the line
%
% Algorithm: <algorithm_name>
%
% then <algorithm_name> is returned as the human-readable algorithm name. 
% Otherwise, if this line cannot be found, the function name is returned.
%
% Author: Tim Mullen 2012, SCCN/INC, UCSD.
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


persistent mvarAlgs;

defaultNameOnly = false;
algName = '';

if nargin==1 && strcmp(varargin{1},'defaultNameOnly')
    defaultNameOnly = true;
end

if nargin==2 && ismember_bc('mfileNameOnly',varargin);
    algName = varargin{2};
    if ~ischar(algName)
        error('SIFT:est_fitMVAR:hlp_getMVARalgorithms:badInput', ...
              'Bad argument pair for ''mfileNameOnly'''); 
    end
end

if ~isempty(mvarAlgs)
    % only do the check once (it's time-consuming)
    % ... retrive "cached" version
    algs = mvarAlgs;
else
    
    algs = {};
        
    % get the names of all mvar_* functions in the /est/mvar folder
    siftroot = hlp_getSiftRoot;
    fpath    = [siftroot filesep 'est' filesep 'mvar' filesep];
    mvarFcns = wildcardsearch(fpath,'*mvar_*.m',true,true);
    mvarFcns = regexprep(mvarFcns,['.*mvar' filesep],'');
    mvarFcns = regexprep(mvarFcns,'\.m','');
    
    % cycle through the list of algorithms
    for k=1:length(mvarFcns)
        
        % initialize variable which will determine whether 
        % we exclude this algorithm (for instance if one or
        % more dependencies are missing)
        skipAlgorithm = false; 
        
        % get the help text (H1) for the algorithm entry function
        preText = hlp_getFcnPreambleText(mvarFcns{k});
        
        % extract the human-readable algorithm name from the text 
        % following the 'Algorithm:' header
        preText = strtrim(preText);
        algName = regexpi(preText,'Algorithm\s*[:]?\s*([^\n]*)','tokens');
        if ~isempty(algName)
            algName = strtrim(algName{1}{1});
        else
            algName = mvarFcns{k};
        end
        
        % search if there are any key dependencies (functions) missing from
        % the path. Key dependencies (m-file names) are listed in the 
        % 'Dependencies:' block of the function HelpText preamble.
        deplist = regexpi(preText,'Dependencies\s*[:]?\s*([^\n]*)','tokens');
        if ~isempty(deplist)
            % parse the list of function dependencies 
            % and extract m-file/function names
            deplist = regexp(deplist{1}{1},'\S*[^ \(\),\s]','match');
            deplist = strtrim(deplist); % just in case

            for dep=1:length(deplist)
                % check if the dependency function exists
                if ~exist(deplist{dep},'file')
                    % ... a critical dependency does not exist
                    % so we have to exclude this algorithm
                    disp_once('WARNING: The MVAR algorithm ''%s'' depends on %s.m, which cannot be located on the path. This algorithm will not be available.',algName,deplist{dep});
                    skipAlgorithm = true;
                    break;
                end
            end
        end
        
        if ~skipAlgorithm
            algs{end+1} = {algName str2func(mvarFcns{k})};
        end
    end
    
    mvarAlgs = algs;
end

% return only the name of the first available algorithm
if defaultNameOnly
    algs = 'Vieira-Morf';
%     algs = algs{1}{1};
elseif ~isempty(algName)
    % get the function handle matching the desired algorithm name
    tmp = cellfun(@(x) fastif(isequal(x{1},algName),x{2},''),algs,'UniformOutput',false);
    tmp(cellfun(@isempty,tmp))=[];
    if isempty(tmp)
        error('SIFT:est_fitMVAR:hlp_getMVARalgorithms:badAlgorithmName', ...
              'Unknown algorithm ''%s''',algName);
    else
        % return function name as a string
        algs = char(tmp{1});
    end
end