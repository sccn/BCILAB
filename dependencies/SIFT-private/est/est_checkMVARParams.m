function [infostring warnstring errstring] = est_checkMVARParams(EEG, params)
%
% This function performs a series of sanity checks on the parameters chosen
% for fitting a VAR model. The results are returned in cell vectors
% indicating information bulletins (notifications), warnings, and/or
% errors.
% 
% Inputs:
%
%   EEG:            EEG data structure
%   params:         a parameter structure as output from est_fitMVAR
%
% Outputs:
%
%   infostring:     cell array of information bulletins (one per cell)
%   warnstring:     cell array of warning bulletins (one per cell)
%   errstring:      cell array of error bulletins (one per cell)
%
% References:
% 
% [1] Korzeniewska, et al (2008). Dynamics of Event-Related Causality in
%     Brain Electrical Activity. Human brain mapping 29:1170?1192 
% [2] Schlogl and Supp, (2006). Analyzing event-related EEG data with 
%     multivariate AR parameters. Prog. in Brain Researc. Neuper and Klimesh, Eds.
% [3] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual.
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% See Also: pop_est_fitMVAR()
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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

srate = EEG.srate;
wlen = srate*params.winlen;  % winlen in samples
T = EEG.CAT.trials; 
p = params.morder;
M = EEG.CAT.nbchan;

% addtional error checks and info
infostring = {};
warnstring = {};
errstring  = {};

% critical error checks
if params.winstep==0
    error('Step size must be greater than 0\n');
end

if floor(params.winlen*1000) > floor((EEG.xmax-EEG.xmin+1)*1000)
    error('Winlen cannot be greater than the trial length (%0.1f sec)\n',EEG.xmax-EEG.xmin);
end


if round(params.winlen*EEG.srate) <= max(params.morder)
    error('The window length must be greater than the model order. Increase your window to at least %0.2f sec\n',(params.morder(end)+1)/EEG.srate);
end


if length(p)==2
    infostring = [infostring sprintf('Two model orders specified [%d %d]\n',p(1),p(2))];
    warnstring = [warnstring sprintf('\tI assume you are providing a [min max] range for model order selection.\n\tI will use p=(%d) for the remaining checks...\n',p(2))];
    errstring  = [errstring {''}];
    p = p(2);
end

% Check if we have a reasonable ratio of datapoints to free parameters
% Optimal ratio derived from [1]
winrat = ((M^2)*p)/(wlen*T);   
infostring = [infostring sprintf('Ratio of number of parameters to datapoints is %0.3f.\n',winrat)];
if winrat > 1
    reclen = ((M^2)*p/T)/srate;
    warnstring = [warnstring sprintf('\tIf using an unconstrained (unregularized) model fitting apprach: The ratio of number of parameters to datapoints must be <= 1.\n\tYour window length must be at least %0.3f sec\n',reclen)];
    errstring  = [errstring {''}];
elseif winrat > 0.1
    reclen = (10*(M^2)*p/T)/srate;
    warnstring = [warnstring sprintf('\tIf using an unconstrained (unregularized) model fitting apprach: For best results, ratio of number of parameters to datapoints should be < 0.1.\n\tI recommend using window length of at least %0.3f sec\n',reclen)];
    errstring  = [errstring {''}];
else
    warnstring = [warnstring {''}];
    errstring  = [errstring {''}];
end

% check if time-frequency principle is violated [2]
TFprod = wlen*sqrt(T);
infostring = [infostring sprintf('Time-Frequency Product is %0.3f. This should be greater than p=%d\n',TFprod,p)];
if wlen<=p/sqrt(T)
    reclen = p/sqrt(T);
    warnstring = [warnstring sprintf('\tTime window does not satisfy the Time-Frequency Uncertainty Principle.\n\tI recommend a window length of at least %0.3f sec\n', reclen)];
    errstring  = [errstring {''}];
else
    warnstring = [warnstring {''}];
    errstring  = [errstring {''}];
end


% notify user of max number of spectral peaks
infostring = [infostring sprintf('Given your model order of p=%d, a maximum of p/2=%0.1f frequency components (spectral peaks) can be estimated for each pair of variables\n',p, p/2)];
warnstring = [warnstring {''}];
errstring  = [errstring {''}];



