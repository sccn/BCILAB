
function [Conn peaks] = hlp_collapseConn(varargin)
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

Conn = arg_extract(varargin,'Conn',2);
if ~isstruct(Conn) || isempty(Conn)
    Conn = struct([]);
end
g = arg_define(varargin, ...
    arg_norep({'Conn','Connectivity'},mandatory,[],'Connectivity structure. Can also be a PConn structure.'), ...
    arg({'connmethods','ConnectivityMethods'},hlp_getConnMethodNames(Conn),hlp_getConnMethodNames(Conn),'Connectivity method names. Cell array of connectivity method names.'), ...
    arg_sub({'coldim','DimensionToCollapse'},{}, ...
    { ...
        arg_subtoggle({'freq','Frequency','Freq'},'off', ...
            { ...
            arg({'range','Range'},[],[],'Value range [min max]. Can also an [N x 2] matrix where each row contains a [min max] range to select.','shape','matrix'), ...
            arg({'method','CollapseMethod'},'net',{'net','mean','peak','peak2d','getrange','maxmag','max','min'},'Collapse method'), ...
            arg({'dim','Dimension'},3,[1 Inf],'Measure dimension. This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to automatically determine dimension from Conn.dims') ...
            arg({'order','Order'},1,{1 2},'Order to apply transformation'), ...
            }, 'Collapse across frequency dimension'), ...
         arg_subtoggle({'time','Time'},'off', ...
            { ...
            arg({'range','Range'},[],[],'Value range [min max]. Can also an [N x 2] matrix where each row contains a [min max] range to select.','shape','matrix') ...
            arg({'method','CollapseMethod'},'net',{'net','mean','peak','peak2d','getrange','maxmag','max','min'},'Collapse method'), ...
            arg({'dim','Dimension'},4,[1 Inf],'Measure dimension. This determines the dimension of Conn.(methodname) to collapse. If empty, we will try to automatically determine dimension from Conn.dims') ...
            arg({'order','Order'},2,{1 2},'Order to apply transformation'), ...
            }, 'Collapse across time dimension'), ...
        },'Dimensions to collapse.'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output') ...
    );
    
    
%     
% 
% % parse inputs
% g = finputcheck(varargin,...
%    {'connmethods', ''          {}          {}; ...      % cell array of connectivity method names. If empty, we use all of them
%     'connThresh'   ''          []          0; ...       % absolute thresholding to apply after significance thresholding. can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
%     'prcThresh'    'real'      [eps 100]   100;...      % top percentile of connections to keep
%     'range'       'real'      []          [];...       % same units as Conn.freqs
%     'range'       'real'      []          [];...       % same units as EEG.CAT.times
%     'sigThresh'    ''          []          [];...       % can be a scalar or a matrix of same size as EEG.CAT.C. Can be logical C(C~=thresh) = 0 or real-valued C(C<thresh)=0
%     'badchans'     'integer'   []          [];...
%     'method'       ''          ''          '';...       % cell array of dimensionality reduction methods to apply in a specified order.
%                                                         % e.g. {'freq','net','time','peak'} will first integrate over freq, then find peak over time.
%                                                         % If method is a string, then this is taken to be the compression global dimension reduction
%                                                         % method applied to all dims in the order {'time','freq',...}
%                                                         % Allowed compression methods are:
%                                                         % 'net'         (integrate)
%                                                         % 'mean'        (average)
%                                                         % 'peak'        (1-D peak search along given dim)
%                                                         % 'peak2d'      (2-D peak search) 
%                                                         % 'getrange'    (extract a data range from given dimension)
%                                                         % 'maxmag'      (1-D max magnitude (max of absval) along given dimension)
%                                                         % 'max'         (1-D max along given dimension)
%                                                         % 'min'         (1-D min along given dimension)
% 
%     'peakWidth'    'integer'   []          2;...
%     'chanlocs'     ''          []          [];...
%     'distRad'      'real'      []          [];...
%     'strictRad'    'boolean'   []          1;...
%     'metric',      'string'    []          'manhattan';...
%     'diags',       'string'    {'off','on'} 'on'; ...
%     'dim',     'real'      []          []; ...
%     'dim',     'real'      []          []; ...
%     'verb',        'boolean'   []          1; ...
%     },'hlp_filterConns','error','quiet');

Conn = g.Conn;

% Determine freq/time dimensions
if ~isempty(g.coldim.time.dim)
    dim = g.coldim.time.dim;
elseif isfield(Conn,'dims')
    dim=find(ismember_bc(Conn.dims,'time'));
else
    error('could not automatically determine the dimension for ''time''. Perhaps it is the 4th dimension?');
end
if ~isempty(g.coldim.freq.dim)
    dim = g.coldim.freq.dim;
elseif isfield(Conn,'dims')
    dim=find(ismember_bc(Conn.dims,'freq'));
else
    error('could not automatically determine the dimension for ''freqs''. Perhaps it is the 3rd dimension?');
end
% Determine default freq and time range
if isempty(g.coldim.freq.range)
    g.coldim.freq.range = [Conn.freqs(1) Conn.freqs(end)];
end
if isempty(g.coldim.time.range)
    g.coldim.time.range = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
end
if length(Conn.freqs)==1
    Conn.freqs = repmat(Conn.freqs,1,2);
end
if length(Conn.erWinCenterTimes)==1
    Conn.erWinCenterTimes = repmat(Conn.erWinCenterTimes,1,2);
    Conn.winCenterTimes = repmat(Conn.winCenterTimes,1,2);
end
if isnan(g.coldim.freq.range)
    g.coldim.freq.range = [];
end
if isnan(g.coldim.time.range)
    g.coldim.time.range = [];
end
frangeidx = getindex(Conn.freqs,g.coldim.freq.range);
trangeidx = getindex(Conn.erWinCenterTimes,g.coldim.time.range);


% determine sequence of collapse operations
colseq  = setdiff_bc(fieldnames(g.coldim),'arg_direct');
colseq  = colseq(cellfun(@(fn) g.coldim.(fn).arg_selection,colseq));
% reorder operations by preferred order
seqorder = zeros(1,length(colseq));
for k=1:length(colseq), seqorder(k) = g.coldim.(colseq{k}).order; end
if length(unique_bc(seqorder))<length(seqorder)
    error('SIFT:hlp_collapseConn', ...
          'Please specify a unique value of ''Order'' for each dimension'); 
end
colseq = hlp_vec(colseq(seqorder));

% dimorder = g.coldim.dimorder;
% dimorder = regexprep(dimorder,{'[fF]requency','[Ff]req'},'freq');
% dimorder = regexprep(dimorder,{'[tT]ime'},'time');

% Collapse each conn method individually
for m=1:length(g.connmethods)
    connmethod = g.connmethods{m};
    
%     if length(size(Conn.(connmethod)))>5
%         error('Connectivity Matrix cannot be greater than 5-D!');
%     end
    
    C = Conn.(connmethod);  
        
    for dd = colseq'
        dim_name = dd{1};
        dim_idx = g.coldim.(dim_name).dim;
        col_method = g.coldim.(dim_name).method;
        
        % collapse time dim
        if (isequal(dim_name,'time')) && ~isempty(g.coldim.time.range)
            if g.verb
                fprintf('applying method=%s to dim=%s\n',dim_name,dim_idx);
            end
            [C peaks.times] = collapseDim(C,'dim',dim_idx,...
                'range',trangeidx,'method',col_method,...
                'dx',Conn.erWinCenterTimes(2)-Conn.erWinCenterTimes(1),'peakWidth',g.peakWidth, ...
                'minpeakheight',-Inf);
        end

        % collapse freq dim
        if ~isempty(g.coldim.freq.range) && (isequal(dim_name,'freq'))
            if g.verb
                fprintf('applying method=%s to dim=%s\n',dim_name,dim_idx);
            end
            [C peaks.freqs] = collapseDim(C,'dim',dim,...
                'range',frangeidx,'method',dim_idx,...
                'dx',Conn.freqs(2)-Conn.freqs(1),'peakWidth',2, ...
                'minpeakheight',-Inf);
        end

    end
            
    
    sz=size(C);
    if all(sz(3:end)==1)
        C = squeeze(C);
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



