function bsplmodel = stat_bsplsm_mkspl(ntimes,nfreqs,timeKnots,freqKnots,verb)
% construct the B-spline basis functions

% ntimes and nfreqs are the number of time and frequency points
% timeKnots is the indices (positions) of knots along the time axis
% freqKnots is the indices (positions) of knots along the freq axis
%
% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
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
% interpolation factor

delta=1/1000; 

if nargin<5
    verb=1;
end

if verb
    fprintf('Constructing b-spline basis functions.\n');
    fprintf('Number of knots:\n');
    fprintf('\tTime: %d\n',length(timeKnots));
    fprintf('\tFreq: %d\n',length(freqKnots));
    fprintf('Working...');
end
   
% Set data parameters
% -------------------------------------------------------------------------
times = 1:ntimes;
freqs = 1:nfreqs;

nTimeKnots  = length(timeKnots)+2;
nFreqKnots  = length(freqKnots)+2;
nTotalKnots = nTimeKnots*nFreqKnots;

time_cont=1+(0:((ntimes-1)/delta))*delta;
freq_cont=1+(0:((nfreqs-1)/delta))*delta;

% Construct the spline basis functions
% -------------------------------------------------------------------------
sp1=cell(nTimeKnots,1);
sp2=cell(nFreqKnots,1);
% construct time-dimension splines
for k=1:nTimeKnots
    sp1{k} = spmak(augknt(timeKnots,4),[zeros(1,k-1) 1 zeros(1,nTimeKnots-k)],nTimeKnots);
end
% construct freq-dimension splines
for k=1:nFreqKnots
    sp2{k} = spmak(augknt(freqKnots,4),[zeros(1,k-1) 1 zeros(1,nFreqKnots-k)],nFreqKnots);
end

% differentiate and evaluate splines
S1=zeros(nTimeKnots);
for h1=1:nTimeKnots
    fnvalh1 = fnval(fnder(sp1{h1},2),time_cont);
    for h2=1:h1
        fnvalh2 = fnval(fnder(sp1{h2},2),time_cont);
        S1(h1,h2) = delta*fnvalh1*fnvalh2';
        S1(h2,h1) = S1(h1,h2);
    end
end

S2=zeros(nFreqKnots);
for h1=1:nFreqKnots
    fnvalh1 = fnval(fnder(sp2{h1},2),freq_cont);
    for h2=1:h1
        fnvalh2 = fnval(fnder(sp2{h2},2),freq_cont);
        S2(h1,h2)=delta*fnvalh1*fnvalh2';
        S2(h2,h1)=S2(h1,h2);
    end
end

% Transform basis into Smooth and Wiggly parts
% -------------------------------------------------------------------------
S1_twdle= kron(S1,eye(nFreqKnots));
S2_twdle= kron(eye(nTimeKnots),S2);
S_star  = S1_twdle+S2_twdle;
U       = eig(S_star);
U_F     = U(:,1:4);
U_R     = U(:,5:nTotalKnots);
S1_star = U_R'*S1_twdle*U_R;
S2_star = U_R'*S2_twdle*U_R;
S1_star = (S1_star'+S1_star)/2;
S2_star = (S2_star'+S2_star)/2;
Phi_F   = zeros(4,ntimes,nfreqs);
Phi_R   = zeros((nTotalKnots-4),ntimes,nfreqs);

for t=1:ntimes
    if verb && ~mod(t,ntimes/10), fprintf('.%d',100*(t/ntimes)); end
    for f=1:nfreqs
        phi_tf=reshape((spcol(augknt(timeKnots,4),4,times(t))' ...
            *(spcol(augknt(freqKnots,4),4,freqs(f))))',nTotalKnots,1);
        Phi_F(:,t,f) = U_F'*phi_tf;
        Phi_R(:,t,f) = U_R'*phi_tf;
    end
end

if verb, fprintf('\n'); end

% Return the spline model struct
% -------------------------------------------------------------------------
bsplmodel = struct('S1_star',   S1_star,    ...
                   'S2_star',   S2_star,    ...
                   'Phi_F',     Phi_F,      ...
                   'Phi_R',     Phi_R,      ...
                   'nTimeKnots',nTimeKnots, ...
                   'nFreqKnots',nFreqKnots, ...
                   'ntimes',    ntimes,     ...
                   'nfreqs',    nfreqs,     ...
                   'class',     'bivariate' ...
                   );
                   

               