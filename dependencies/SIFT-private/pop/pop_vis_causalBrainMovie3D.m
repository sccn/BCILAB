
function [ALLEEG cfg handles] = pop_vis_causalBrainMovie3D(ALLEEG,typeproc,varargin)
%
% Create an interactive 3D BrainMovie from a connectivity matrix. See [1]
% for more details on the Interactive BrainMovie3D. 
% This function generates an interactive GUI. Default GUI options can be
% pre-specified as <name, value> pairs as listed in vis_causalBrainMovie3D().
% Alternately, an arg configuration structure can be provided in lieu of
% name, value pairs (see Example below).
%
% Input: 
%   
%   ALLEEG:     Array of EEGLAB structures containing Connectivity matrix
%
% Optional:
% 
%   Name, Value pairs as defined in vis_causalBrainMovie3D()
%
% Output:
%
%   cfg:    an argument configuration structure which can be used to
%           repopulate GUI fields or in a command-line call to
%           vis_causalBrainMovie3D(). 
%
% Example:  >> cfg = pop_vis_causalBrainMovie3D(ALLEEG);    % now select some options from GUI and click "Make Movie!"
%           >> vis_causalBrainMovie3D(ALLEEG,cfg);          % Movie will now be recreated from command-line 
%                                                           % exactly as specified in GUI 
%
% See Also: vis_causalBrainMovie3D()
% 
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen and Arnaud Delorme, 2010, SCCN/INC, UCSD. 
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

handles = [];

if length(varargin) == 1 && isempty(varargin{1})
    varargin = {};
end

if nargin<2
    typeproc = 0;
end

fcnName     = strrep(mfilename,'pop_','');
fcnHandle   = str2func(fcnName);

% check the dataset
res = hlp_checkeegset(ALLEEG,{'conn'});
if ~isempty(res)
    error(['SIFT:' fcnName],res{1});
end

if any(cellfun(@isempty,{ALLEEG.dipfit}))
    error('SIFT:vis_causalBrainMovie3D:NeedDipfit', ...
        'BrainMovie3D requires source or channel locations stored in EEG.dipfit');
end

if isfield(ALLEEG(1).CAT.configs,fcnName)
    % get default configuration (from prior use) and merge with varargin
    varargin = [hlp_struct2varargin(ALLEEG(1).CAT.configs.(fcnName)) varargin];
end


for cnd=1:length(ALLEEG);
    % make node labels a row array
    ALLEEG(cnd).CAT.curComponentNames = ALLEEG(cnd).CAT.curComponentNames(:)';
end

if strcmpi(typeproc,'nogui')
    % get the config from function
    cfg = arg_tovals(arg_report('rich',fcnHandle,[{'EEG',ALLEEG(1),'Conn',ALLEEG(1).CAT.Conn},varargin]),false);
else
    % render the GUI
    [PGh figh] = feval(['gui_' fcnName],ALLEEG(1),ALLEEG(1).CAT.Conn,varargin{:});
    
    if isempty(PGh)
        % user chose to cancel
        cfg = [];
        return;
    end
    
    % get the specification of the PropertyGrid
    ps = PGh.GetPropertySpecification;
    cfg = arg_tovals(ps,false);
end

drawnow;

if strcmpi(typeproc,'cfg_only')
    return;
end

% execute the low-level function
for cnd=1:length(ALLEEG)
    handles = feval(fcnHandle,'EEG',ALLEEG(cnd),'Conn',ALLEEG(1).CAT.Conn,cfg);
    
    if ~isempty(cfg)
        % store the configuration structure
        ALLEEG(cnd).CAT.configs.(fcnName) = cfg;
    end
end
