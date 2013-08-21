function [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
% Helper function to convert structure params to variables used by the
% various routines - also performs checks to ensure that parameters are
% defined; returns default values if they are not defined.
%
% Usage: [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
%
% Inputs:
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%           - optional
%             tapers : precalculated tapers from dpss or in the one of the following
%                       forms:  
%                       (1) A numeric vector [TW K] where TW is the
%                           time-bandwidth product and K is the number of
%                           tapers to be used (less than or equal to
%                           2TW-1). 
%                       (2) A numeric vector [W T p] where W is the
%                           bandwidth, T is the duration of the data and p 
%                           is an integer such that 2TW-p tapers are used. In
%                           this form there is no default i.e. to specify
%                           the bandwidth, you have to specify T and p as
%                           well. Note that the units of W and T have to be
%			                consistent: if W is in Hz, T must be in seconds
% 			                and vice versa. Note that these units must also
%			                be consistent with the units of params.Fs: W can
%		    	            be in Hz if and only if params.Fs is in Hz.
%                           The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Outputs: 
% The fields listed above as well as the struct params. The fields are used
% by some routines and the struct is used by others. Though returning both
% involves overhead, it is a safer, simpler thing to do.

if ~isfield(params,'tapers') || isempty(params.tapers);  %If the tapers don't exist
     display('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers=[3 5];
end;
if ~isempty(params) && length(params.tapers)==3 
    % Compute timebandwidth product
    TW = params.tapers(2)*params.tapers(1);
    % Compute number of tapers
    K  = floor(2*TW - params.tapers(3));
    params.tapers = [TW  K];
end

if ~isfield(params,'pad') || isempty(params.pad);
    params.pad=0;
end;
if ~isfield(params,'Fs') || isempty(params.Fs);
    params.Fs=1;
end;
if ~isfield(params,'fpass') || isempty(params.fpass);
    params.fpass=[0 params.Fs/2];
end;
if ~isfield(params,'err') || isempty(params.err);
    params.err=0;
end;
if ~isfield(params,'trialave') || isempty(params.trialave);
    params.trialave=0;
end;

tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;
err=params.err;
trialave=params.trialave;