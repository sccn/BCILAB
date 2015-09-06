function hlpText = hlp_buildSimHelpText()
% return a formatted help text string which contains help text definitions
% for each of the simulation functions
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
% first get the list of function names

sims = hlp_getSimExamples;

hlpText = '';

% for each algorithm...
for k=1:length(sims)
    % ... get the human-readable name for the algorithm
    simFullName = sims{k}{1};
    % ... get the function m-file name ...
    simFcnName = hlp_getSimExamples('mfileNameOnly',simFullName);
    % ... get the preamble for the function ...
    try preamble = hlp_getFcnPreambleText(simFcnName);
    catch err
        disp(err.message);
        preamble = 'No Description Available';
    end
    % ... and insert into the hlpText
    hlpText = sprintf([hlpText '\n\n' ...
                       '---------------------------------------\n' ...
                       simFullName '\n' ...
                       '---------------------------------------\n' ...
                       preamble
                       ]);
end