% pop_readcog() - import data from a CogniScan data file.
%                   
% Usage:
%   >> EEGOUT = pop_readcog( EEG ); % pop-up a data entry window 
%   >> EEGOUT = pop_readcog( 'key', val,...); % no pop-up window
%
% Graphic interface:
%   "EEGLAB dataset name" - [Edit box] name for the new dataset. 
%                  Command line equivalent: 'setname'
%   "Data (.bin) file" - [Edit box] Data file to import to EEGLAB.
%                  Command line equivalent: 'data'
%   "Channels" - [Edit box] Channels to read 
%                  Command line equivalent: 'channel'
%                  0 -> Defaults to number of channels in file
%   "Gain Correction" - [Edit box] Whether to apply gain to data
%                       Command line equivalent: 'gaincorrect'
%                       0 -> No; 1 -> Yes (default)
%   "Time points per epoch" - [Edit box] Number of points per data epoch.
%                  Irrelevant for continuous data. Command line equivalent: 'pnts'
%                  0 -> all data
%   "Optional epoch start time" - [Edit box] Command line equivalent: 'xmin'
%   "List of epoch onsets" - [Edit box] Onset time(s) (in samples) of
%                            epoch(s)
%                            0 -> beginning of file
%   "Channel locations file or array" - [Edit box] see readlocs() help for
%                  data channel format. Command line equivalent: 'chanlocs'
%   "ICA weights array or text file" - [edit box] Use this option to import
%                  ICA weights from other decompositions (for instance: same
%                  data, different conditions). To use the ICA weights from
%                  another loaded dataset (n) enter 'ALLEEG(n).icaweights'
%                  in this edit box. Command line equivalent: 'icaweights'
%   "ICA sphere array or text file" - [edit box] Import an ICA sphering matrix. 
%                  For computational reasons, an ICA decomposition may be defined 
%                  by a sphere matrix and an unmixing (weight) matrix (above).
%                  To use the ICA weights from another loaded dataset (n)
%                  enter "ALLEEG(n).icasphere". If no sphering matrix, enter 
%                  'eye(EEG.nbchan)'. Command line equivalent: 'icasphere'.
% Optional inputs:
%   'setname'    - ['string'] Name of the new EEGLAB dataset
%   'data'       - ['varname'|'filename'] Data variable or file name to import.
%   'chanlocs'   - ['varname'|'filename'] Import a file containing electrode 
%                  locations (see >> help readlocs for file format).
%   'channel'    - Channels to read; default all channels; 'nbchan'=length('channel')
%   'xmin'       - Starting time in seconds.
%   'offset'     - Starting time of epochs in samples
%   'pnts'       - Number of data frames (time points) per data epoch (epoched data only).
%   'icaweights' - ICA weight matrix. 
%   'icasphere'  - ICA sphering matrix (if [], eye(nchans)).
% 
% Outputs:
%   EEGOUT      - modified EEG dataset structure
%
% Adapted from pop_importdata
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson and Lucas Parra
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

function [EEGOUT, com] = pop_readcog(varargin);

com = '';
EEGOUT = eeg_emptyset;
if nargin < 1                 % if several arguments, assign values 
   % popup window parameters	
   % -----------------------
   geometry    = { [2 0.1 0.8 0.5] [1.3 0.8 .8 0.5] [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] ...
           [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] ...
           [2 0.1 0.8 0.5] [1.5 0.6 0.8 0.5] [2 0.1 0.8 0.5] [2 0.1 0.8 0.5] };
   commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', tagtest), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
	commandsetfiletype = [ 'filename = get( findobj(''parent'', gcbf, ''tag'', ''globfile''), ''string'');' ...
					'tmpext = findstr(filename,''.'');' ...
					'tmpext = lower(filename(tmpext(end)+1:end));' ...
					'switch tmpext, ' ...
					'  case ''mat'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',5);' ...
					'  case ''fdt'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',3);' ...
					'  case ''txt'', set(findobj(gcbf,''tag'', ''loclist''), ''value'',2);' ...
					'end; clear tmpext filename;' ];
    uilist = { ...
         { 'Style', 'text', 'string', 'EEGLAB dataset name (optional):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', '' }, { }...
         ...
         { 'Style', 'text', 'string', 'Data (.bin) file', 'horizontalalignment', 'right', 'fontweight', 'bold' }, ...
         { }, { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'globfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', ...
		   [ 'tagtest = ''globfile'';' commandload commandsetfiletype ] }, ...
         ...
         { 'Style', 'text', 'string', 'Channels (0->set from data):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, {},  { 'Style', 'edit', 'string', '0' }, { } ...
         { 'Style', 'text', 'string', 'Gain Correction (0->No; 1->Yes):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, {},  { 'Style', 'edit', 'string', '1' }, { } ...
         { 'Style', 'text', 'string', 'Time points per epoch (0=continuous data):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { },  { 'Style', 'edit', 'string', num2str(EEGOUT.pnts) }, { } ...
         { 'Style', 'text', 'string', 'Optional epoch start time for data epochs (sec):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold' }, { }, { 'Style', 'edit', 'string', num2str(EEGOUT.xmin) }, { },...
         { 'Style', 'text', 'string', 'List of epoch onsets (samples; 0=continuous):', 'horizontalalignment', 'right', ...
		   'fontweight', 'bold'}, { }, { 'Style', 'edit', 'string', '0' }, { },...
         ...
         { 'Style', 'text', 'string', 'Channel locations file or array:', 'horizontalalignment', 'right', 'fontweight', ...
		   'bold' }, {'Style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''readlocs'');' }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'chanfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''chanfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA weights array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'weightfile' }, ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''weightfile'';' commandload ] }, ...
         ...
         { 'Style', 'text', 'string', 'ICA sphere array or text file (if any):', 'horizontalalignment', 'right' }, { }, ...
         { 'Style', 'edit', 'string', '', 'horizontalalignment', 'left', 'tag',  'sphfile' } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', [ 'tagtest = ''sphfile'';' commandload ] } };

    results = inputgui( geometry, uilist, 'pophelp(''pop_readcog'');', 'Import Cogniscan dataset info -- pop_readcog()');
    if length(results) == 0, return; end;

	args = {};
	if ~isempty( results{1} ), args = { args{:}, 'setname', results{1} }; end;
	if ~isempty( results{2} ), args = { args{:}, 'dataformat', 'bin' }; end;
	if ~isempty( results{2} ) , args = { args{:}, 'data', results{2} }; end;
	if ~isempty( results{3} ) , args = { args{:}, 'channels', str2num(results{3}) }; end;
    if ~isempty( results{4} ) , args = { args{:}, 'gaincorrect', str2num(results{4}) }; end;
	if ~isempty( results{5} ) , args = { args{:}, 'pnts', str2num(results{5}) }; end;
    if ~isempty( results{6} ) , args = { args{:}, 'xmin', str2num(results{6}) };
        if abs(str2num(results{6})) > 10,
            fprintf('WARNING: are you sure the epoch start time (%3.2f) is in seconds\n');
        end;
    end;
    i=7;
    if ~isempty( results{i} ) , args = { args{:}, 'offset', results{i} }; end;
	if ~isempty( results{i+1} ) , args = { args{:}, 'chanlocs' , results{i+1} }; end;
	if ~isempty( results{i+2} ),  args = { args{:}, 'icaweights', results{i+2} }; end;
	if ~isempty( results{i+3} ) , args = { args{:}, 'icasphere', results{i+3} }; end;

    
    
    
else % no interactive inputs
    EEGOUT = pop_editset(EEGOUT, varargin{:});
    args = varargin;
    for index=1:2:length(args)
        if ~isempty(inputname(index+1)) & ~isstr(args{index+1}) & length(args{index+1})>1, 
			args{index+1} = inputname(index+1); 
        end;
    end;                
end;

% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('Setevent: wrong syntax in function arguments'); return; end;
else
    g = [];
end;

g.offset=['[' g.offset ']'];
g.offset = evalin('base', g.offset, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', g.offset)' );

% assigning values
% ----------------
tmpfields = fieldnames(g);
for curfield = tmpfields'
    switch lower(curfield{1})
        case 'offset' , ;
        case {'dataformat' }, ; % do nothing now
        case 'setname' , EEGOUT.setname = getfield(g, {1}, curfield{1});
        case 'pnts'    , EEGOUT.pnts = getfield(g, {1}, curfield{1});
        case 'comments', EEGOUT.comments = getfield(g, {1}, curfield{1});
	    case 'averef'  , disp('The ''averef'' argument is obsolete; use function pop_reref() instead');
        case 'channels'  , EEGOUT.nbchan = length(getfield(g, {1}, curfield{1}));

%       case 'srate'   , EEGOUT.srate = getfield(g, {1}, curfield{1});
        case 'chanlocs', varname = getfield(g, {1}, curfield{1});
                         if isempty(varname)
                             EEGOUT.chanlocs = [];
                         elseif exist( varname ) == 2
                            fprintf('Pop_editset: channel locations file ''%s'' found\n', varname); 
                            EEGOUT.chanlocs = readlocs(varname);
                         else
                            EEGOUT.chanlocs = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
                         end;
        case 'icaweights', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: ICA weight matrix file ''%s'' found\n', varname); 
							try, EEGOUT.icaweights = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch, fprintf('Pop_editset warning: error while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icaweights = [];
							 else
								 EEGOUT.icaweights = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
							 end;
						 end;
                         if ~isempty(EEGOUT.icaweights) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
        case 'icasphere', varname = getfield(g, {1}, curfield{1});
                         if exist( varname ) == 2
                            fprintf('Pop_editset: ICA sphere matrix file ''%s'' found\n', varname); 
                            try, EEGOUT.icasphere = load(varname, '-ascii');
								EEGOUT.icawinv = [];
                            catch, fprintf('Pop_editset warning: erro while reading filename ''%s'' for ICA weight matrix\n', varname); 
                            end;
                         else
							 if isempty(varname) 
								 EEGOUT.icasphere = [];
							 else
								 EEGOUT.icasphere = evalin('base', varname, 'fprintf(''Pop_editset warning: variable ''''%s'''' not found, ignoring\n'', varname)' );
								 EEGOUT.icawinv = [];
							 end;
                         end;
                         if ~isempty(EEGOUT.icaweights) & isempty(EEGOUT.icasphere)
                            EEGOUT.icasphere = eye(size(EEGOUT.icaweights,2));
                         end;
	    case 'data'    , varname = getfield(g, {1}, curfield{1});
                         if ~isempty(g.channels), 
                             EEGOUT.nbchan=length(g.channels); 
                             EEGOUT.comments = pop_comments({EEGOUT.comments 'Channels: ' num2str(g.channels)},[],[],1);

                         end 
                         if EEGOUT.nbchan == 0,
								EEGOUT.nbchan=[]; g.channels=[];
                                fprintf(['Channels unspecified, reading all channels...\n']);
						 end;   
                         if g.pnts==0,
                             g.pnts==[];
                         end
                         if g.offset==0,
                             g.offset=[];
                         end
                         
                         try,
                             [EEGOUT.data,D,N,EEGOUT.srate,gain] = readcogdata(varname(1:end-4), g.channels, g.gaincorrect, g.pnts, g.offset);
                            
                             EEGOUT.comments = pop_comments({EEGOUT.comments 'Gain: ' num2str(gain)},[],[],1);
                             EEGOUT.pnts=g.pnts;
                             if (isempty(EEGOUT.pnts)|(EEGOUT.pnts==0)),
                                 EEGOUT.pnts=N;
                             end
                             if (isempty(g.channels)|(g.channels==0)),
                                 EEGOUT.nbchan=D; g.channels=1:D;
                             end
                             
                             EEGOUT.trials=max(1,length(g.offset));
                             EEGOUT.data=1e6.*reshape(EEGOUT.data,[EEGOUT.nbchan EEGOUT.pnts EEGOUT.trials]);
                            
                             
                          catch, error(['Pop_readcog error: cannot read Cogniscan data file ''' varname ''' ']); 
                          end;
        case 'xmin'    , if exist('EEG.xmin'), oldxmin = EEG.xmin; end
                         EEGOUT.xmin = getfield(g, {1}, curfield{1});
                         EEGOUT.xmax = (EEGOUT.pnts-1)/EEGOUT.srate+EEGOUT.xmin;
                         if exist('EEG.event'),
                            if ~isempty(EEG.event)
                             if nargin < 2
                                if ~popask( ['Warning: changing the starting point of epochs will' 10 'lead to recomputing epoch event latencies ?'] )
                                    error('Pop_editset: transformation cancelled by user');
                                end;
                             end;
                             if isfield(EEG.event, 'latency')
                                for index = 1:length(EEG.event)
                                    EEG.event(index).latency = EEG.event(index).latency - (EEG.xmin-oldxmin)*EEG.srate;
                                end;
                             end;       
                            end;
                        end;
         case 'gaincorrect' , ;
             
         otherwise, error(['Pop_editset error: unrecognized field ''' curfield{1} '''']); 
    end;
end;

    % generate the output command
    % ---------------------------
    com = '';
    for i=1:2:length(args)
        if ~isempty( args{i+1} )
            if isstr( args{i+1} ) com = sprintf('%s, ''%s'', ''%s''', com, args{i}, char(args{i+1}) );
            else                  com = sprintf('%s, ''%s'', [%s]', com, args{i}, num2str(args{i+1}) );
            end;
        else
            com = sprintf('%s, ''%s'', []', com, args{i} );
        end;
    end;


    com = [ 'EEG = pop_readcog(' com(2: end) ');']; 
    

return;
