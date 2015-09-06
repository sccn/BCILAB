function [C, Q] = procCosp(DATA, WINDOW_TYPE, WINDOW_SIZE, OVERLAP, FREQUENCY_RANGE, PHASE_CORRECTION)
% [C, Q] = procCosp(DATA, WINDOW_TYPE, WINDOW_SIZE, OVERLAP, FREQUENCY_RANGE)
%
% This function computes cospectra and quadrature spectra on multivariate 
% data using Welch's method. For a better estimation, each tappering window 
% can be phase corrected according to its length and its position in data. 
%
% Inputs:
% - DATA            --> matrix with the data (M samples by N channels matrix)
% - WINDOW_TYPE     --> (optionnal) string with the type of the window
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples 
% - OVERLAP         --> (optionnal) percentage of the overlapping window (from 0 to 0.99)
% - FREQUENCY_RANGE --> (optionnal) range of interesting frequencies to keep, 
%                                   format must be [sample_rate min_frec max_frec]
% - PHASE_CORRECTION--> (optionnal) correction of phase shift intrinsic to Welch's windows. 
%                                   (1: correction, 0: no correction)
%
% Outputs:
% - C        --> cospectral matrix [N x N x f]
% - Q        --> quadrature matrix [N x N x f]
%
% History
% First version: Jonas Chatel-Goldman @ GIPSA-Lab, 01/12/2011


% default UI display
if(nargin < 9)  UIdisplay  = 0;  end
% default correction
if(nargin < 6)  PHASE_CORRECTION = 0; end
% default frequency range
if(nargin < 5)  FREQUENCY_RANGE = []; end
% default overlapping
if(nargin < 4)  OVERLAP = .75;  end
% default windows size
if(nargin < 3)  WINDOW_SIZE = 128;  end
% default windows type
if(nargin < 2)  WINDOW_TYPE = 'hanning';  end

% Adjust the window to the next pow of 2 (for the FFT)
WINDOW_SIZE = 2 ^ nextpow2(WINDOW_SIZE);
% Calculate the number of frequencies
number_freqs = (WINDOW_SIZE / 2); % +1;

s_data = size(DATA);


% UI DISPLAY
if(UIdisplay)  h= waitbar(0,'Processing cospectrum...');  end

%% Create the window from the inputs and preallocates memory
switch lower(WINDOW_TYPE)
    case 'hamming'
        win = hamming(WINDOW_SIZE);
    case 'hann'
        win = hann(WINDOW_SIZE);
    case 'hanning'
        win = hanning(WINDOW_SIZE);
    case 'blackman'
        win = blackman(WINDOW_SIZE);
    case 'barthannwin'
        win = barthannwin(WINDOW_SIZE);
    case 'blackmanharris'
        win = blackmanharris(WINDOW_SIZE);
    case 'bohmanwin'
        win = bohmanwin(WINDOW_SIZE);
    case 'chebwin'
        win = chebwin(WINDOW_SIZE);
    case 'gausswin'
        win = gausswin(WINDOW_SIZE);
    case 'kaiser'
        win = kaiser(WINDOW_SIZE);
    case 'nuttallwin'
        win = nuttallwin(WINDOW_SIZE);
    case 'parzenwin'
        win = parzenwin(WINDOW_SIZE);
    case 'tukeywin'
        win = tukeywin(WINDOW_SIZE);
    case 'rectangular'
        win = ones(WINDOW_SIZE,1);
    case 'flattopwin'
        win = flattopwin(WINDOW_SIZE);
    case 'welch'
         %win = welchwin(WINDOW_SIZE, 0);
        Nd2=(WINDOW_SIZE-1)/2;
        for i=1 : WINDOW_SIZE
            win(i)=1-((i-Nd2)/Nd2)^2;
        end
        win = win';
    otherwise
        disp('WARNING - Unknown window type!');
        disp('Switching to default window type : Hamming Window');
        win = hamming(WINDOW_SIZE);
end

if s_data(1) < WINDOW_SIZE
    error('Cospectra computation: data segment too short for this window length, consider changing frequency resolution');
end


% Calculate the number of windows because of the overlapping
if(OVERLAP == 0)
    number_windows = floor(s_data(1) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data(1)/WINDOW_SIZE); 
    number_windows = 1 + (nbFullWin-1)/(1-OVERLAP) + floor((s_data(1)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end


% pre-allocation of memory 
S = zeros(s_data(2), s_data(2), number_freqs);

%% Loop on all frequencies
for window_ix = 1 : number_windows
    
    % UI display
    if(UIdisplay)  waitbar(window_ix/number_windows,h); end
    
    % time markers to select the data
    t1 = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                           % marker of the end of the time window
    
    % select current window and apodize it   
    cdata = DATA(t1:t2, :) .* (win*ones(1,s_data(2)));

    % FFT calculation
    fdata = fft(cdata,WINDOW_SIZE) ./ number_windows; % complex data
    if(PHASE_CORRECTION == 1)
        fdata = fdata.*(exp(-sqrt(-1)*t1*( 0:(WINDOW_SIZE-1) )'/WINDOW_SIZE*2*pi)*ones(1,s_data(2)));
    end
    
    % Complexe cospectrum averaging over windows
    for f = 1 : number_freqs
        S(:,:,f) = S(:,:,f) + (fdata(f,:)' * fdata(f,:));
    end
end

C = real(S);
Q = imag(S);


% Adjust Frequency range to specified range (in case it is a parameter)
if ~isempty(FREQUENCY_RANGE)     
    if  FREQUENCY_RANGE(3) <= (FREQUENCY_RANGE(1) / 2)        
        Faxis = [0:(1/number_freqs):1-(1/number_freqs)];
        FindexMin = find(Faxis >= FREQUENCY_RANGE(2)/(FREQUENCY_RANGE(1)/2), 1);
        FindexMax = find(Faxis >= FREQUENCY_RANGE(3)/(FREQUENCY_RANGE(1)/2), 1);
        C = C(:, :, FindexMin:FindexMax);
        Q = Q(:, :, FindexMin:FindexMax);

    else         
        % Assertion on frequency, to comply with Nyquist Theorem
        error('Frequency range exceeds data sampling rate in consideration of Shannon''s limitation.')
    end
end

% UI display
if(UIdisplay) close(h); end