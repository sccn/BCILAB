
function [Conn peaks] = hlp_filterConns(Conn, varargin)
%
% Apply a set of selection rules to a connectivity matrix
%
% Inputs:
%
%   Conn        Connectivity structure with subfields Conn.(connmethod)
%               containing [nchs x nchs x num_freqs x num_times]
%               connectivity matrices
%
% Optional:
%
%       'connThresh'    [real]. Absolute thresholding to apply after
%                               significance thresholding. Can be be a 
%                               scalar or a matrix of same size as 
%                               EEG.CAT.Conn.(connmethod). Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
%       'prcThresh'     -     [real] The upper percent [1-100] of filtered
%                             connections to keep {def: 100}
%       'method'              cell array of dimensionality reduction methods to apply in a specified order. e.g.,
%                             {'freq','net','time','peak'} will first integrate over freq, then find peak over time.
%                             If method is a string, then this is taken to be the compression global dimension reduction
%                             method applied to all dims in the order {'time','freq',...}.
%                             allowed methods:  'mean', 'net', 'peak', 'getrange'
%        ...
%
    
    
    
% Outputs:
%
%   Conn        filtered connectivity structure
%   peaks       frequency peaks (if method = 'peak')
%
% See Also: est_mvarConnectivity()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
%
%
% Author: Tim Mullen 2008, 2010, SCCN/INC, UCSD.
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


% ----------------------------------------------------------
% TODO:
% need to provide ability to specify order in which dim-compression
% methods are applied. e.g., {3,'net',4,'peak'} would first
% integrate over the third dimension (frequency) and then find the peaks
% over the 4th dim (time)
% ONE OPTION:  allow 'method' = a cell array as above
% ----------------------------------------------------------
% Key point, need to have EEG.CAT.dims = {'chans','chans','freq','times'}
% or some such...
%


if nargin<1
    help hlp_filterConns;
    return;
end


% parse inputs
g = finputcheck(varargin,...
   {'connmethods', ''          {}          {}; ...      % cell array of connectivity method names. If empty, we use all of them
    'connThresh'   ''          []          0; ...       % absolute thresholding to apply after significance thresholding. can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
    'prcThresh'    'real'      [eps 100]   100;...      % top percentile of connections to keep
    'frange'       'real'      []          [];...       % same units as Conn.freqs
    'trange'       'real'      []          [];...       % same units as EEG.CAT.times
    'sigThresh'    ''          []          [];...       % can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
    'badchans'     'integer'   []          [];...
    'method'       ''          ''          '';...       % cell array of dimensionality reduction methods to apply in a specified order.
                                                        % e.g. {'freq','net','time','peak'} will first integrate over freq, then find peak over time.
                                                        % If method is a string, then this is taken to be the compression global dimension reduction
                                                        % method applied to all dims in the order {'time','freq',...}
                                                        % Allowed compression methods are:
                                                        % 'net'         (integrate)
                                                        % 'mean'        (average)
                                                        % 'peak'        (1-D peak search along given dim)
                                                        % 'peak2d'      (2-D peak search) 
                                                        % 'getrange'    (extract a data range from given dimension)
                                                        % 'maxmag'      (1-D max magnitude (max of absval) along given dimension)
                                                        % 'max'         (1-D max along given dimension)
                                                        % 'min'         (1-D min along given dimension)

    'peakWidth'    'integer'   []          2;...
    'chanlocs'     ''          []          [];...
    'distRad'      'real'      []          [];...
    'strictRad'    'boolean'   []          1;...
    'metric',      'string'    []          'manhattan';...
    'diags',       'string'    {'off','on'} 'on'; ...
    'freqdim',     'real'      []          []; ...
    'timedim',     'real'      []          []; ...
    'verb',        'boolean'   []          1; ...
    },'hlp_filterConns','error','quiet');

        
if ischar(g)
    error(g);
end

if isempty(g.connmethods)
    g.connmethods = hlp_getConnMethodNames(Conn);
end

nmethods = 3;

if ~isempty(g.timedim)
    timedim = g.timedim;
elseif isfield(Conn,'dims')
    timedim=find(ismember_bc(Conn.dims,'time'));
else
    timedim = [];
end

if ~isempty(g.freqdim)
    freqdim = g.freqdim;
elseif isfield(Conn,'dims')
    freqdim=find(ismember_bc(Conn.dims,'freq'));
else
    freqdim = [];
end

if isempty(freqdim)
    if g.verb
        disp('warning: no ''freq'' entry found in Conn.dims -- assuming dim 3');
    end
    freqdim = 3;
end
if isempty(timedim)
    if g.verb
        disp('warning: no ''time'' entry found in Conn.dims -- assuming dim 4');
    end
    timedim = 4;
end

if ischar(g.method),
    ONEMETHOD = true;
    g.method = {'time',g.method,'freq',g.method};
else
    ONEMETHOD = false;
end

if isempty(g.frange)
    g.frange = [Conn.freqs(1) Conn.freqs(end)];
end

if isempty(g.trange)
    g.trange = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
end

if length(Conn.freqs)==1
    Conn.freqs = repmat(Conn.freqs,1,2);
end

if length(Conn.erWinCenterTimes)==1
    Conn.erWinCenterTimes = repmat(Conn.erWinCenterTimes,1,2);
    Conn.winCenterTimes = repmat(Conn.winCenterTimes,1,2);
end

if isnan(g.frange)
    g.frange = [];
end

if isnan(g.trange)
    g.trange = [];
end

frangeidx = getindex(Conn.freqs,g.frange);
trangeidx = getindex(Conn.erWinCenterTimes,g.trange);

for m=1:length(g.connmethods)
    connmethod = g.connmethods{m};
    
    
    if length(size(Conn.(connmethod)))>5
        error('Connectivity Matrix cannot be greater than 4-D!');
    end
    
    C = Conn.(connmethod);  % make copy of connectivity matrix
    
    
    % apply significance thresholding
    if ~isempty(g.sigThresh)
        if islogical(g.sigThresh) && isequal(size(C),size(g.sigThresh))
            C(~g.sigThresh)=0;
        else
            C(C<g.sigThresh) = 0;
        end
    end
    
    if ~all(cellfun(@isempty,g.method(2:2:end))) && ~strcmpi(g.method{2},'thresh')
        if ~isempty(g.trange) && ~isempty(g.frange) && ONEMETHOD && isequal(g.method{2},'peak2d')
            % apply 2D peak identification
            if g.verb
                fprintf('applying method=%s to dim=%s and %s\n',g.method{2},g.method{1},g.method{3});
            end
            C = collapseDim(C,'dim',4,'range',trangeidx, ...
                'method','getrange');     % extract time range
            C = collapseDim(C,'dim',3,'range',frangeidx, ...
                'method','getrange');     % extract freq range
            [C peaks] = collapseDim(C,'dim',4,'range',[], 'method','peak2d', ...
                'peakWidth',g.peakWidth,'minpeakheight',g.connThresh); % find peaks
        else
            
            for mm=1:2:length(g.method)
                
                curmethod = g.method{mm};
                
                if ~ismember_bc(curmethod,{'time','freq'})
                    error('unknown method %s',curmethod);
                end
                
                % collapse time dim
                if ~isempty(g.trange) && (isequal(curmethod,'time'))
                    if g.verb
                        fprintf('applying method=%s to dim=%s\n',g.method{mm+1},g.method{mm});
                    end
                    [C peaks.times] = collapseDim(C,'dim',timedim,...
                        'range',trangeidx,'method',g.method{mm+1},...
                        'dx',Conn.erWinCenterTimes(2)-Conn.erWinCenterTimes(1),'peakWidth',g.peakWidth, ...
                        'minpeakheight',g.connThresh);
                end
                
                % collapse freq dim
                if ~isempty(g.frange) && (isequal(curmethod,'freq'))
                    if g.verb
                        fprintf('applying method=%s to dim=%s\n',g.method{mm+1},g.method{mm});
                    end
                    [C peaks.freqs] = collapseDim(C,'dim',freqdim,...
                        'range',frangeidx,'method',g.method{mm+1},...
                        'dx',Conn.freqs(2)-Conn.freqs(1),'peakWidth',g.peakWidth, ...
                        'minpeakheight',g.connThresh);
                end
                
            end
            
        end
    end
    
    sz=size(C);
    if all(sz(3:end)==1)
        C = squeeze(C);
    end
    
    % remove diagonals
    if strcmpi(g.diags,'off')
        for q=1:size(C,3)
            C(:,:,q)=C(:,:,q).*~eye(size(C(:,:,q)));
        end
    end
    
    % apply connectivity threshold
    if g.connThresh
        C(C<g.connThresh) = 0;
    end
    
    % apply percentile threshold
    if ~isempty(g.prcThresh) && g.prcThresh<100
        % get specified percentile
        prc = prctile(C(:),100-g.prcThresh);
        C(C<prc)=0;
    end
    
    Conn.(connmethod) = C;
    
end

if ~isempty(frangeidx) && size(C,3)==1
    % frequency dimension has been collapsed
    % replace with median frequency within collapsed interval
    Conn.collapsedFreqs = Conn.freqs(frangeidx(1):frangeidx(end));
    Conn.freqs          = median(Conn.collapsedFreqs);
else
    Conn.freqs          = Conn.freqs(frangeidx(1):frangeidx(end));
end

if ~isempty(trangeidx) &&size(C,4)==1
    % time dimension has been collapsed
    % replace with median time point within collapsed interval
    Conn.collapsedTimes     = Conn.erWinCenterTimes(trangeidx(1):trangeidx(end));
    Conn.erWinCenterTimes   = median(Conn.collapsedTimes);
    Conn.winCenterTimes     = median(Conn.winCenterTimes(trangeidx(1):trangeidx(end)));
else
    Conn.erWinCenterTimes   = Conn.erWinCenterTimes(trangeidx(1):trangeidx(end));
    Conn.winCenterTimes     = Conn.winCenterTimes(trangeidx(1):trangeidx(end));
end



