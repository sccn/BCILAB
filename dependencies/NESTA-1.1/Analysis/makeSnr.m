function sigma = makeSnr( x, snr, quant )
% sigma = makeSnr( x, snrLevel)
%   finds the variance "sigma" appropriate for white-noise
%   so that x vs x + sigma*randn has SNR of snrLevel
%  snrLevel = Inf is acceptable (sigma will be zero)
%
% sigma = makeSnr( x, snrLevel, quantization)
%   where quantization is true, attempts to mimick
%   the effects of a finite-precision ADC
%   i.e. instead of a function of power of signal,
%   it's now a function of peak-to-peak amplitude
%   (e.g. l_infinity norm)
%   If quantization if false, then reverts to usual SNR

rms = @(x) norm(x)/sqrt(length(x));

if nargin < 3, quant = false; end

if snr == Inf
    sigma = 0;
else     
    if quant
        sigma = 10^( log10(norm(x,'inf')) - 1/20*snr );
    else
        sigma = 10^( log10(rms(x)) - 1/20*snr );
    end
end
