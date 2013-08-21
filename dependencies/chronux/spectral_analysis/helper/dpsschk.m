function [tapers,eigs]=dpsschk(tapers,N,Fs)
% Helper function to calculate tapers and, if precalculated tapers are supplied, 
% to check that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
% Usage: tapers=dpsschk(tapers,N,Fs)
% Inputs:
% tapers        (tapers in the form of: 
%                                   (i) precalculated tapers or,
%                                   (ii) [NW K] - time-bandwidth product, number of tapers) 
%
% N             (number of samples)
% Fs            (sampling frequency - this is required for nomalization of
%                                     tapers: we need tapers to be such
%                                     that integral of the square of each taper equals 1
%                                     dpss computes tapers such that the
%                                     SUM of squares equals 1 - so we need
%                                     to multiply the dpss computed tapers
%                                     by sqrt(Fs) to get the right
%                                     normalization)
% Outputs: 
% tapers        (calculated or precalculated tapers)
% eigs          (eigenvalues) 

if nargin < 3; error('Need all arguments'); end

% wrapper for caching of results in RAM
if nargout == 1
    tapers = hlp_microcache('multitaper_dpsschk',@do_compute,tapers,N,Fs);
else
    [tapers,eigs] = hlp_microcache('multitaper_dpsschk',@do_compute,tapers,N,Fs);
end

% do the actual computation (cached on disk for fast lookup)
function [tapers,eigs] = do_compute(tapers,N,Fs)
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=hlp_diskcache('filterdesign',@dpss,N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
