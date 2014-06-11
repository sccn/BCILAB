% pop_loadxdf() - Load an XDF file (*.xdf or *.xdfz).
%                 (pop out window if no arguments)
%
% Usage:
%   >> [EEG] = pop_loadxdf;
%   >> [EEG] = pop_loadxdf( filename, 'key', 'val', ...);
%
% Graphic interface:
%
%   "Stream name to import" - [edit box] specify name of stream to import; if nonempty, only the 
%                             stream with the given name will be imported (otherwise the stream type
%                             will be used to determine what stream to import)
%                             Command line equivalent in eeg_load_xdf: 'streamname'
%   "Stream type to import" - [edit box] specify content type of stream to import
%                             see http://code.google.com/p/xdf/wiki/MetaData (bottom) for content types
%                             Command line equivalent in eeg_load_xdf: 'streamtype'
%   "Exclude marker stream(s)" - [edit box] specify names of marker streams to skip; this is in 
%                                MATLAB cell array syntax, e.g. {'MyVideoMarkers','SyncStream001'}
%                                Command line equivalent in eeg_load_xdf: 'exclude_markerstreams'
%
% Inputs:
%   filename                   - file name
%
% Optional inputs:
%   'streamname'               - name of stream to import (if omitted, streamtype takes precedence)
%   'streamtype'               - type of stream to import (default: 'EEG')
%   'exclude_markerstreams'    - cell array of marker stream names that should be excluded from import
%   Same as eeg_load_xdf() function.
%
% Outputs:
%   [EEG]                       - EEGLAB data structure
%
% Note:
% This script is based on pop_loadcnt.m to make it compatible and easy to use in
% EEGLab.
%
% Author: Christian Kothe, Swartz Center for Computational Neuroscience, UCSD, 2012
%
% See also: eeglab(), eeg_load_xdf(), load_xdf()
%

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2012 Christian Kothe, ckothe@cusd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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


function [EEG, command]=pop_loadxdf(filename, varargin);

command = '';
filepath = '';
EEG=[];

if nargin < 1

	% ask user
	[filename, filepath] = uigetfile('*.xdf;*.xdfz', 'Choose an XDF file -- pop_loadxdf()');
    drawnow;
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
    uigeom     = { [1 0.5] [1 0.5] [1 0.5] 0.13};
    uilist   = { { 'style' 'text' 'string' 'Stream name to import:' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'string' 'Stream type to import:' } ...
                 { 'style' 'edit' 'string' 'EEG' } ...
                 { 'style' 'text' 'string' 'Exclude marker streams(s):' } ...
                 { 'style' 'edit' 'string' '{}' } {}};

	result = inputgui(uigeom, uilist, 'pophelp(''pop_loadxdf'')', 'Load an XDF file');
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    options = [];
    if ~isempty(result{1}),
        options = [options ', ''streamname'', ''' result{1} '''']; end
    if ~isempty(result{2}),
        options = [options ', ''streamtype'', ''' result{2} '''']; end
    if ~isempty(result{3}),        
        options = [options ', ''exclude_markerstreams'', ' result{3} '']; end
else
	options = vararg2str(varargin);
end;

% load data
% ----------
if exist('filepath','var')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;

fprintf('Now importing...');
if nargin > 0    
    EEG = eeg_load_xdf(fullFileName, varargin{:});
else
	eval( [ 'EEG = eeg_load_xdf( fullFileName ' options ');' ]);
end;
fprintf('done.\n');

EEG = eeg_checkset(EEG);

if length(options) > 2
    command = sprintf('EEG = pop_loadxdf(''%s'' %s);',fullFileName, options);
else
    command = sprintf('EEG = pop_loadxdf(''%s'');',fullFileName);
end;

