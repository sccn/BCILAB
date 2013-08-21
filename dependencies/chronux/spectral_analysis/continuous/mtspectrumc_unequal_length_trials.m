function [ S, f, Serr ]= mtspectrumc_unequal_length_trials( data, movingwin, params, sMarkers )

% This routine computes the multi-taper spectrum for a given set of unequal length segments. It is
% based on modifications to the Chronux routines. The segments are continuously structured in the 
% data matrix, with the segment boundaries given by markers. Below,
% movingwin is used in a non-overlaping way to partition each segment into
% various windows. Th spectrum is evaluated for each window, and then the
% window spectrum estimates averaged. Further averaging is conducted by
% repeating the process for each segment. 
%
% Inputs: 
%
%   data = data( samples, channels )- here segments must be stacked
%   as explained in the email 
%   movingwin = [window winstep] i.e length of moving
%              window and step size. Note that units here have
%              to be consistent with units of Fs. If Fs=1 (ie normalized)
%              then [window winstep]should be in samples, or else if Fs is
%              unnormalized then they should be in time (secs). 
%   sMarkers = N x 2 array of segment start & stop marks. sMarkers(n, 1) = start
%           sample index; sMarkers(n,2) = stop sample index for the nth segment
%   params = see Chronux help on mtspecgramc
%
% Output:
%
%       S       frequency x channels
%       f       frequencies x 1
%       Serr    (error bars) only for err(1)>=1
%
%

iwAvg = 1; % 0=no weighted average, 1=weighted average
debug = 1; % will display intermediate calcs. 

if nargin < 2; error('avgSpectrum:: Need data and window parameters'); end;
if nargin < 3; params=[]; end;
if isempty( sMarkers ), error( 'avgSpectrum:: Need Markers...' ); end
[ tapers, pad, Fs, fpass, err, trialave, params ] = getparams( params );
if nargout > 2 && err(1)==0; 
%   Cannot compute error bars with err(1)=0. change params and run again.
    error('avgSpectrum:: When Serr is desired, err(1) has to be non-zero.');
end;

% Set moving window parameters to no-overlapping
if abs(movingwin(2) - movingwin(1)) >= 1e-6, disp( 'avgSpectrum:: Warming: Window parameters for averaging should be non-overlapping. Set movingwin(2) = movingwin(1).' ); end

wLength = round( Fs * movingwin(1) ); % number of samples in window
wStep = round( movingwin(2) * Fs ); % number of samples to step through

% Check whether window lengths satify segment length > NW/2
if ( wLength < 2*tapers(1) ), error( 'avgSpectrum:: movingwin(1) > 2*tapers(1)' ); end

% Left align segment markers for easier coding
sM = ones( size( sMarkers, 1 ), 2 ); 
sM( :, 2 ) = sMarkers( :, 2 ) - sMarkers( :, 1 ) + 1;

% min-max segments 
Nmax = max( sM(:,2) ); Nmin = min( sM(:,2) );
if ( Nmin < 2*tapers(1) ), error( 'avgSpectrum:: Smallest segment length > 2*tapers(1). Change taper settings' ); end

% max time-sample length will be the window length. 
nfft = 2^( nextpow2( wLength ) + pad );
[ f, findx ] = getfgrid( Fs, nfft, fpass); 

% Precompute all the tapers
sTapers = tapers;
sTapers = dpsschk( sTapers, wLength, Fs ); % compute tapers for window length

nChannels = size( data, 2 ); 
nSegments = size( sMarkers, 1 );

if debug
    disp( ['Window Length = ' num2str(wLength)] );
    disp( ['Window Step = ' num2str(wStep)] );
    disp( ' ' );
end

s = zeros( length(f), nChannels );
serr = zeros( 2, length(f), nChannels );
S = zeros( length(f), nChannels );
Serr = zeros( 2, length(f), nChannels );
nWins = 0;
for sg = 1 : nSegments
    % Window lengths & steps fixed above
    % For the given segment, compute the positions & number of windows
    N = sM(sg,2); 
    wStartPos = 1 : wStep : ( N - wLength + 1 );
    nWindows = length( wStartPos );
    if nWindows
        nWins = nWins + nWindows; % for averaging purposes

        w=zeros(nWindows,2);
        for n = 1 : nWindows
            w(n,:) = [ wStartPos(n), (wStartPos(n) + wLength - 1) ]; % nWindows x 2. just like segment end points
        end

        % Shift window limits back to original sample-stamps
        w(:, 1) = w(:,1) + (sMarkers( sg, 1 ) - 1);
        w(:, 2) = w(:,2) + (sMarkers( sg, 1 ) - 1);

        if debug
            disp( ['Segment Start/Stop = ' num2str( w(1,1) ) ' ' num2str( w(end,2) ) ] );
            disp( ['Min / Max Window Positions = ' num2str( min(w(:,1)) ) ' ' num2str( max(w(:,1)) ) ] );
            disp( ['Total Number of Windows = ' num2str(nWindows) ]);
            disp( ' ' );
        end

        % Pile up window segments similar to segment pileup
        wData = zeros( wLength, nChannels, nWindows ); %initialize to avoid fragmentation
        for n = 1:nWindows
            %wData( :, :, n ) = detrend( data( w(n,1):w(n,2), : ), 'constant' );
            wData( :, :, n ) = detrend( data( w(n,1):w(n,2), : ) );
        end

        % J1 = frequency x taper x nWindows
        % J2 = frequency x taper x nWindows x nChannels
        J2 = zeros( length(f), tapers(2), nWindows, nChannels ); J2 = complex( J2, J2 );
        for c = 1 : nChannels
            J1 = mtfftc( squeeze(wData( :, c, : )), sTapers, nfft, Fs ); % FFT for the tapered data
            J2( :, :, :, c ) = J1(findx,:,:);
        end
        % J2 = frequency x taper x nWindows x nChannels
        % Inner mean = Average over tapers => frequency x nWindows x nChannels
        % Outer mean = Average over windows => frequency x nChannels
        dim1 = [length(f), nWindows, nChannels];
        dim2 = [length(f), nChannels];
        % s = frequency x nChannels
        s = reshape( squeeze( mean( reshape( squeeze( mean( conj(J2).*J2, 2 ) ), dim1), 2 ) ), dim2 );

        % Now treat the various "windowed data" as "trials"
        % serr = 2 x frequency x channels. Output from specerr = 2 x frequency x 1
        for c = 1 : nChannels
            serr( :, :, c ) = specerr( squeeze( s(:, c ) ), squeeze( J2(:,:,:, c ) ), err, 1 );
        end
        
        if iwAvg
            % Segment Weighted error estimates.
            S = S + nWindows*s;
            Serr = Serr + nWindows*serr;
        else
            S = S + s;
            Serr = Serr + serr;
        end

    else
        if debug, disp(['avgSpectrum:: Zero windows for segment: ' num2str(sg) ]); end
    end
end

% Segment Weighted error estimates.
% Only over those that had non-zero windows
if nWins && iwAvg
    S=S/nWins; Serr=Serr/nWins;
end
if ~nWins
    if debug, disp(['avgCoherence:: No segment long enough with movingwin parameters found. Reduce movingwin.' ]); end
end




