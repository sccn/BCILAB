function rgbdata = cmplx2rgb(cmplxdata, logmag)
%CMPLX2RGB         Convert complex matrix to an imageable RGB form.
%   RGBDATA = CMPLX2RGB(CMPLXDATA) converts the [M x N] CMPLXDATA to an
%   [M x N x 3] truecolor matrix in RGB color coordinates.  The mapping
%   goes through HSB space, where the hue is determined by the data phase,
%   the brightness is determined by the data magnitude, and the saturation
%   is 0.6.  The RGB values are scaled to be between 0 and 1 (set to 1 if
%   all data points have the same magnitude).
%
%   RGBDATA = CMPLX2RGB(CMPLXDATA, 'log') uses the log of the data
%   magnitude to determine brightness.  Zero values are set to NaN.

% Argument checking.
if (nargin < 2)
    logmag = 0;
elseif (strcmp(logmag, 'log'))
	logmag = 1;
else
	error('Unknown argument.');
end

% Compute amplitude and phase.
amplitude = abs(cmplxdata);
phase = angle(cmplxdata) + pi;   % shift from [-Pi, Pi] to [0, 2*Pi].
if (logmag)
    old = warning('off', 'MATLAB:log:logOfZero');
    amplitude = log10(amplitude);
    amplitude(isinf(amplitude)) = NaN;
    warning(old);
end
ampmin = min(amplitude(:));
ampmax = max(amplitude(:));

% Construct HSV color representation.
hsvdata = zeros([size(cmplxdata) 1]);
hsvdata(:,:,1) = phase / (2 * pi);
hsvdata(:,:,2) = 0.6;
if (abs(ampmax-ampmin) > 2*eps)   % must have variation above roundoff error
    hsvdata(:,:,3) = (amplitude - ampmin) ./ (ampmax - ampmin);
else
    hsvdata(:,:,3) = 1;
end

% Finally, convert HSV to RGB.
rgbdata = hsv2rgb(hsvdata);