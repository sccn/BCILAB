
% collapseDim() - collapse a range of values across a chosen dimension
%                 using the chosen method (e.g., average, peak, integral).
%
% Inputs:
%   C           -   [Nchs x Nchs x DIM3 x Dim4 x ... ] connectivity matrix.
%
% Optional:
%   'dim'         -   the dimension to collapse {def: 3}
%   'range'       -   [min-idx max-idx] index range to collapse. 
%                     If range is a single value, than the layer for that  
%                     index only is returned. {def: collapse entire dim}
%   'method'      -   ['mean' | 'net' | 'peak']. 'net' applies the trapz() 
%                     integral ('dx' option must be set) {def: 'mean'}
%                 -   if method=@fcn, then function fcn() will be applied.
%                     This function must take a matrix as first argument
%                     and a dimension across which to evaluate as the
%                     second arg.
%   'dx'          -   [real] Required if 'method'='net'. This is the 
%                     spacing increment (in same units as dataset)for trapz
%                     e.g., if integrating across frequency with freq.
%                     spacing of 0.5Hz, then dx=0.5;
%   'peakWidth'   -   [real] The width (in points) over which to identify
%                     a local peak. e.g., if peakWidth is 10, C(n,n,f) is
%                     only a peak if C(n,n,f:f-5)<C(n,n,f)>C(n,n,f:f+5)
%   'numPeaks'    -   [real] maximum number of peaks to return [if
%                     method='peaks' or 'peaks2d']
% Output:
%   C       -   collapsed [Nchs x Nchs x ... x 1 x ...] matrix with 
%               dimension 'dim' reduced to a singleton
%   peaks   -   [Nchs x Nchs x ... x (numdims(C)-1) x Npeaks] 
%               matrix containing indicies of peaks along the collapsed 
%               dimension.
%
% Example:  collapse C(:,:,10:50) using integral with unit spacing
%           C2 = collapseDim(C,'dim',3,'range',[10 50],'method','net','dx',1)
%
% SEE ALSO
% localmax(), CAT
%
% Copyright (c) 2008 by Tim Mullen
function [C peaks] = collapseDim(C,varargin)
    
if nargin<2
    help collapseFreqs;
    return;
end

peaks = [];

range = [];
method = 'mean';
peakWidth = 20;
dim = 3;
dx = 0;
mph = -Inf;
numPeaks = 1;
keepedges = 1;
extent = 0;

for i=1:2:length(varargin)-1
    switch varargin{i}
        case 'range',       range = varargin{i+1};
        case 'method',      method = varargin{i+1};
        case 'dx',          dx = varargin{i+1};
        case 'dim',         dim = varargin{i+1};
        case 'peakWidth',   peakWidth = varargin{i+1};
        case 'numPeaks',    numPeaks = varargin{i+1};
        case 'minpeakheight', mph = varargin{i+1};
        case 'keepedges',   keepedges = varargin{i+1};
        case 'extent',      extent = varargin{i+1};
    end
end

if size(C,dim)==1
    return;
end

maxrange = size(C,dim);

if dim>length(size(C))
    error('The matrix does not have this many dimensions!');
end

% collapse over full range?
if isempty(range)
    range = [1 maxrange];
end

if length(range)>2
    error('range must be singleton or a 2-element [min max] vector');
end

% handle when single freq. input
if length(range)==1,
    C = getRange(C,dim,[range range]);
    return;
end

% range checking
if range(1)<1 || range(2)>maxrange
    error(sprintf(['range out of bounds\n' ...
                   'range must be between %d and %d'],...
                   1,maxrange));
   return;
end


switch lower(method)
    case 'mean'
        C = nan_mean(getRange(C,dim,range),dim);
    case 'net'
        if ~dx, error('''dx'' must be specified for ''net'' method'); end
        C = trapz(getRange(C,dim,range),dim);
        C = C.*dx;
    case 'peak'     % 1D peak detection
        [C peaks] = findConnPeaks(getRange(C,dim,range),dim, mph,peakWidth,numPeaks,keepedges,extent);
    case 'peak2d'   % 2D peak detection
%         C = getRange(C,dim,range);
        newC = zeros([size(C,1) size(C,2)],'single');
        peaks.freqs = zeros(size(C,1),size(C,2),numPeaks);
        peaks.times = zeros(size(C,1),size(C,2),numPeaks);
        if issymmetric(C)    % only estimate peaks for upper triangle (faster)
            for i=1:size(C,1)
                for j=i:size(C,2)
                    [pkval peaks.freqs(i,j,:) peaks.times(i,j,:)] = findpeaksND(squeeze(C(i,j,:,:)),peakWidth, true,numPeaks);   %max([size(C,3) size(C,4)])
                    newC(i,j) = pkval; %max(pkval);  % only return the maximum peak
                    newC(j,i) = newC(i,j);   % symmetric, so fill in lower triangle too
                end
            end
        else
            for i=1:size(C,1)
                for j=1:size(C,2)  % all channels

                    % new stuff
                    [pk] = imregionalmax(squeeze(C(i,j,:,:)));
                    if all(pk(:))
                        continue;  % no peak here
                    else
                        cc = squeeze(C(i,j,:,:));
                        [rows cols] = size(cc);
                        if extent
                            [ii jj] = find(pk);
                            for pp=1:length(ii)  % for each peak
                                % check if it's extent-neighborhood is
                                % greater than mph
                                
                                if all(cc(setdiff_bc(max(1,ii-extent):min(rows,ii+extent),ii),setdiff_bc(max(1,jj-extent):min(cols,jj+extent),jj))<=mph)
                                    pk(ii,jj)=0;
                                end
                            end
                            if ~any(pk(:)), continue; end
                        end
                        
                        cc = cc.*pk;
                        [newC(i,j) ii] = max(cc(:));
                        [peaks.freqs(i,j) peaks.times(i,j)] = ind2sub(size(cc),ii);
                        
                    end
                    
                end
            end
        end
        C = newC;
        
    case 'maxmag'
        % 1-D max magnitude (max of absval) along dimension dim
        [C peaks] = max(abs(getRange(C,dim,range)),[],dim);
    case 'max'
        % 1-D max along dimension dim
        [C peaks] = max(getRange(C,dim,range),[],dim);
    case 'min'
        % 1-D min along dimension dim
        [C peaks] = min(getRange(C,dim,range),[],dim);
    case {'getrange' 'shrinkonly'}
        C = getRange(C,dim,range);
    otherwise
        try
            C = feval(method,getRange(C,dim,range),dim);
        catch
            error('could not evaluate function');
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract a range of values from a specified dim
function C = getRange(C,dim,range)

numdims = length(size(C));

st = '';
for i=1:numdims
    st = fastif(dim==i,[st ' range(1):range(2),'],[st ':,']);
end
st(end) = [];

eval(['C = C( ' st ');']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Cout pidxs] = findConnPeaks(C,dim,mph,mpd,npeaks,keepedges,extent)
    % C is [nchs nchs dim3 dim4]
    % dim is the dimension to find 1-d peaks along (this dimension will be
    % collapsed to a singleton)
    % output is same size as C with singleton dimension dim 
    % mph is the minimum peak heights
    % mpd is the minimum distance between peaks
    % npeaks is the max number of peaks to return
    % keepedges determines whether or not to count a max at the series edge
    %   as a peak
    % extent determines the number of pixels adjecent to a peak that must
    % be greater than mph
    
    if nargin<6
        keepedges = 0;
    end
    
    if nargin<7
        extent = 0;
    end
    
    warning off
    
    verb = 0;
    
    if nargin<5
        npeaks = 10;    % number of top peaks to store in pidxs
    end
    
    
    % permute C so we are finding peaks over last (4th) dimension
    permutation = [1 2 fastif(dim==3,[4 3],[3 4])];
    C=permute(C,permutation);
    
    [rows cols nd3 nd4] = size(C);
    pidxs = zeros(rows,cols,nd3,npeaks);
    Cout = zeros(rows,cols,nd3);
    
    if nd3>10
        fprintf('warning: searching for peaks over %d series -- this may take a while...\n',nd3);
        verb=1;
    end
    if verb, h=waitbar(0,'searching for peaks...'); cur=1; end
    for kk=1:nd3
        for ch1=1:rows
            for ch2=1:cols
                series = [squeeze(C(ch1,ch2,kk,:))];
                if keepedges==0
                    [peaks locs] = findpeaks(series,length(series),extent,0,mph);
                else
                    [peaks locs] = findpeaks([0 ; series ; 0],length(series),extent,0,mph);
                end
                locs = round(locs);
                
                if ~isempty(locs)
                    pidxs(ch1,ch2,kk,1:min(length(locs),npeaks)) = ...          % keep indices of npeaks largest peak 
                        locs(1:min(length(locs),npeaks))-keepedges;             % -1 to counteract zero-padding
                    Cout(ch1,ch2,kk)=peaks(1);                                  % store largest peak value
                else
                    pidxs(ch1,ch2,kk,:)=nan;
                    Cout(ch1,ch2,kk)=nan;
                end
                
                if verb, waitbar(cur/(nd3*rows*cols),h); cur=cur+1; end
                
            end
        end
        
        if ~mod(kk,10), fprintf('.'); end
    end
    
    if verb, close(h); end
    warning on;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sym = issymmetric(C)
% check if C is symmetric along the first two dimensions

sym = 1;

sz = size(C);

for i=1:sz(1)
    for j=1:sz(2)
        if ~isequal(C(i,j,:,:,:,:),C(j,i,:,:,:,:))
            sym = 0;
            return;
        end
    end
end


