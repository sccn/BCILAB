
function [data g] = pre_selectchans(varargin)
% 
% Normalization is an important operation for improving stationarity when
% fitting models to ensembles of data (data with multiple
% realizations/trials). This function performs ensemble or temporal 
% normalization. Ensemble normalization removes the ensemble mean 
% (e.g. ERP) and divides by the ensemble standard deviation. Temporal 
% normalization removes the mean from each trial and divides each trial 
% by its standard deviation. 
%
%
% Inputs:
%
%   data:        data structure [channels x time x trials]
%
% Optional:     <'Name',Value> pairs
%
%     VerbosityLevel: Verbosity level. 0 = no output, 1 = text, 2 = graphical  
%                     Possible values: 0,1,2                                   
%                     Default value  : 0                                       
%                     Input Data Type: real number (double)                    
% 
%     Method:         Normalize windows across time, ensemble, or both         
%                     Possible values: 'ensemble','time'                       
%                     Default value  : 'time','ensemble'                       
%                     Input Data Type: boolean     
% Outputs:
%
%   EEG:        processed EEG structure
%   g:          argument specification structure
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
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



g = arg_define(varargin, ...
    arg_norep({'signal','Data'}), ...
    arg({'channels','Channels'}, [], [], 'Cell array of channel names to retain.','type','cellstr','shape','row'));

subset = set_chanid(signal,channels);
if ~isequal(subset,1:signal.nbchan)
    data = pop_select(signal,'channel',subset,'sorttrial','off'); end %#ok<*NODEF>


g = arg_define([0 1], varargin, ...
        arg_norep('data',mandatory), ...
        arg({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
        arg({'method','Method'},{'time','ensemble'},{'time','ensemble'},'Normalize windows across time, ensemble, or both','type','logical'));
    
% commit data variable to workspace
data = g.data;
g = rmfield(g,'data');

[nchs pnts ntr] = size(data);

if ischar(g.method)
    g.method = {g.method};
end
       
if g.verb==2
    multiWaitbar('Normalizing','Reset','Color',hlp_getNextUniqueColor);
end

for k=1:length(g.method)
    if g.verb==2
        multiWaitbar('Normalizing',k/length(g.method));
    end
    switch lower(g.method{k})
        case 'ensemble'
            % pointwise subtract ensemble mean (over trials)
            % divide by ensemble stdev
            if ntr==1
                if g.verb, fprintf('Multiple trials not available, ignoring ensemble normalization\n'); end
            else
                if g.verb, fprintf('Normalizing data across ensemble...\n'); end
                normdata = bsxfun(@minus,data,mean(data,3));
                normdata = bsxfun(@rdivide,normdata,std(data,0,3));
            end
        case 'time'
            % pointwise subtract trial mean
            % divide by trial stdev
            if g.verb, fprintf('Normalizing data across time...\n'); end
            normdata = bsxfun(@minus,data,mean(data,2));
            normdata = bsxfun(@rdivide,normdata,std(data,0,2));
    end
    
    data = normdata;
end

if g.verb==2
    multiWaitbar('Normalizing','Close');
end