function [Cxy, F] = cmtm(X,Y,NW,Fs)
%CMTM              Coherence function estimate using the multitaper method.
%   Cxy = CMTM(X,Y) estimates the coherence between two equal length
%   vectors vectors X and Y using Thomson's multitaper method.  The
%   coherence is a complex function of frequency, with magnitude between 0
%   and 1, that estimates
%            E[X*(f) Y(f)] / sqrt(E[X*(f) X(f)] E[Y*(f) Y(f)])
%   where * denotes complex conjugation.  **Note that this is different
%   from the matlab function COHERE, which returns the magnitude squared
%   of this value.**
%
%   Cxy = CMTM(X,Y,NW) specifies the time-bandwidth product for the
%   discrete prolate spheroidal sequences (DPSS) is specified by NW; this
%   value also determines the number of tapers as (2*NW-1).  If not
%   given, the default NW is 4.
%   
%   Cxy = CMTM(X,Y,NW,Fs) specifies a sampling frequency Fs.
%
%   [Cxy,F] = CMTM(...) also returns the vector of frequencies at which
%   the coherence is estimated.  If Fs is specified, this vector ranges
%   from [0,Fs/2] and the units are the same as Fs.  If Fs is not
%   specified, F ranges from [0,1/2] and the units are relative to the 
%   sampling frequency.
%
%   CMTM(...) without output arguments plots the magnitude-squared
%   and phase of the coherence (in two subplots) in the current figure.

%%%%% Debugging test code:
% Fs = 100;
% t = 0:1/Fs:(32-1/Fs);
% X = sin(2*pi*5*t);  Y = sin(2*pi*5*t) + 0.2*randn(size(X));
% cmtm(X,Y,4,Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 2), error('Three arguments are required.');  end
if (nargin < 3), NW = 4;  end;
if (nargin < 4), Fs = 1;  end;

if (size(X,1) > 1), X = X';  end;
if (size(Y,1) > 1), Y = Y';  end;
[Mx,Nx] = size(X);    [My,Ny] = size(Y);
if (Mx~=1 || My~=1),  error('Matrix arguments are not supported.');  end;
if (Nx ~= Ny),        error('A Coherence estimate requires two vectors of equal length.');  end;
if (length(NW) ~= 1), error('NW must be a scalar.');  end;
if (length(Fs) ~= 1), error('Fs must be a scalar.');  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make Tapers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 2*NW-1;
[E,V] = dpss(Nx,NW,P);
E = E';

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Coherence %%%%%%%%%%%%%%%%%%%%%%%
midway = floor(Nx/2)+1;
Xf = fft(E .* repmat(X,P,1),[],2);      Xf = Xf(:,1:midway);     Pxx = conj(Xf).*Xf;    
Yf = fft(E .* repmat(Y,P,1),[],2);      Yf = Yf(:,1:midway);     Pyy = conj(Yf).*Yf;
Cxy = (conj(Xf) .* Yf) ./ sqrt((Pxx.*Pyy));
Cxy = mean(Cxy,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot if Needed %%%%%%%%%%%%%%%%%%%%%%%%%%
F = (0:(1/Nx):0.5)*Fs;

if (nargout == 0)
    subplot(2,1,1);  plot(F,abs(Cxy));    ylabel('Magnitude Squared Coherence');  
    subplot(2,1,2);  plot(F,angle(Cxy));  ylabel('Coherence Phase');   xlabel('Frequency');

    clear Cxy
end

if (nargout < 2)
    clear F;
end;
