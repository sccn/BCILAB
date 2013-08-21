function [outdata,outstate] = asr_process(data,srate,state,windowlen,lookahead,stepsize,maxdims,maxmem,usegpu)
% Processing function for the Artifact Subspace Reconstruction (ASR) method.
% [Data,State] = asr_process(Data,SamplingRate,State,WindowLength,LookAhead,StepSize,MaxDimensions,MaxMemory,UseGPU)
%
% This function is used to clean multi-channel signal using the ASR method. The required inputs are 
% the data matrix, the sampling rate of the data, and the filter state (as initialized by
% asr_calibrate). If the data is used on successive chunks of data, the output state of the previous 
% call to asr_process should be passed in.
%
% In:
%   Data : Chunk of data to process [#channels x #samples]. This is a chunk of data, assumed to be
%          a continuation of the data that was passed in during the last call to asr_process (if
%          any). The data should be *zero-mean* (e.g., high-pass filtered the same way as for
%          asr_calibrate).
%   
%   SamplingRate : sampling rate of the data in Hz (e.g., 250.0)
%
%   State : initial filter state (determined by asr_calibrate or from previous call to asr_process)
%
%   WindowLength : Length of the statistcs window, in seconds (e.g., 0.5). This should not be much
%                  longer than the time scale over which artifacts persist, but the number of samples 
%                  in the window should not be smaller than 1.5x the number of channels. Default: 0.5
%
%   LookAhead : Amount of look-ahead that the algorithm should use. Since the processing is causal,
%               the output signal will be delayed by this amount. This value is in seconds and should
%               be between 0 (no lookahead) and WindowLength/2 (optimal lookahead). The recommended
%               value is WindowLength/2. Default: WindowLength/2
%
%   StepSize : The statistics will be updated every this many samples. The larger this is, the faster 
%              the algorithm will be. The value must not be larger than WindowLength*SamplingRate.
%              The minimum value is 1 (update for every sample) while a good value is 1/3 of a second.
%              Note that an update is always performed also on the first and last sample of the data
%              chunk. Default: 32
%
%   MaxDimensions : Maximum dimensionality of artifacts to remove. Up to this many dimensions (or up 
%                   to this fraction of dimensions) can be removed for a given data segment. If the
%                   algorithm needs to tolerate extreme artifacts a higher value than the default
%                   may be used (the maximum fraction is 1.0). Default 0.66
%
%   MaxMemory : The maximum amount of memory used by the algorithm when processing a long chunk with
%               many channels, in MB. The recommended value is at least 256. To run on the GPU, use
%               the amount of memory available to your GPU here (needs the parallel computing toolbox).
%               default: min(5000,1/2 * free memory in MB). Using smaller amounts of memory leads to
%               longer running times.
%
%   UseGPU : Whether to run on the GPU. This makes sense for offline processing if you have a a card
%            with enough memory and good double-precision performance (e.g., NVIDIA GTX Titan or
%            K20). Note that for this to work you need to have the Parallel Computing toolbox.
%            Default: false
%
% Out:
%   Data : cleaned data chunk (same length as input but delayed by LookAhead samples)
%
%   State : final filter state (can be passed in for subsequent calls)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-31

% UC Copyright Notice
% This software is Copyright (C) 2013 The Regents of the University of California. All Rights Reserved.
% 
% Permission to copy, modify, and distribute this software and its documentation for educational,
% research and non-profit purposes, without fee, and without a written agreement is hereby granted,
% provided that the above copyright notice, this paragraph and the following three paragraphs appear
% in all copies.
% 
% Permission to make commercial use of this software may be obtained by contacting:
% Technology Transfer Office
% 9500 Gilman Drive, Mail Code 0910
% University of California
% La Jolla, CA 92093-0910
% (858) 534-5815
% invent@ucsd.edu 
% 
% This software program and documentation are copyrighted by The Regents of the University of
% California. The software program and documentation are supplied "as is", without any accompanying
% services from The Regents. The Regents does not warrant that the operation of the program will be
% uninterrupted or error-free. The end-user understands that the program was developed for research
% purposes and is advised not to rely exclusively on the program for any reason.
% 
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
% CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS.

if nargin < 4 || isempty(windowlen) 
    windowlen = 0.5; end
windowlen = max(windowlen,1.5*size(data,1)/srate);
if nargin < 5 || isempty(lookahead)
    lookahead = windowlen/2; end
if nargin < 6 || isempty(stepsize)
    stepsize = 32; end
if nargin < 7 || isempty(maxdims)
    maxdims = 0.66; end
if nargin < 9 || isempty(usegpu)
    usegpu = false; end
if nargin < 8 || isempty(maxmem)
    if usegpu
        dev = gpuDevice(); maxmem = dev.FreeMemory/2^20; 
    else
        maxmem = hlp_memfree/(2^21);
    end
end
if maxdims < 1
    maxdims = round(size(data,1)*maxdims); end
if isempty(data)
    outdata = data; outstate = state; return; end

[C,S] = size(data);
N = round(windowlen*srate);
P = round(lookahead*srate);
[T,M,A,B] = deal(state.T,state.M,state.A,state.B);

% initialize prior filter state by extrapolating available data into the past (if necessary)
if isempty(state.carry)
    state.carry = repmat(2*data(:,1),1,P) - data(:,1+mod(((P+1):-1:2)-1,S)); end

data = [state.carry data];
data(~isfinite(data(:))) = 0;

% split up the total sample range into k chunks that will fit in memory
splits = ceil((C*C*S*8*8 + C*C*8*S/stepsize + C*S*8*2 + S*8*5) / (maxmem*1024*1024 - C*C*P*8*3));
if splits > 1
    fprintf('Now cleaning data in %i blocks',splits); end

for i=1:splits
    range = 1+floor((i-1)*S/splits) : min(S,floor(i*S/splits));
    if ~isempty(range)
        % get spectrally shaped data X for statistics computation (range shifted by lookahead)
        [X,state.iir] = filter(B,A,double(data(:,range+P)),state.iir,2);
        % move it to the GPU if applicable
        if usegpu && length(range) > 1000
            try X = gpuArray(X); catch,end; end
        % compute running mean covariance (assuming a zero-mean signal)
        [Xcov,state.cov] = moving_average(N,reshape(bsxfun(@times,reshape(X,1,C,[]),reshape(X,C,1,[])),C*C,[]),state.cov);
        % extract the subset of time points at which we intend to update
        update_at = min(stepsize:stepsize:(size(Xcov,2)+stepsize-1),size(Xcov,2));
        % if there is no previous R (from the end of the last chunk), we estimate it right at the first sample
        if isempty(state.last_R)
            update_at = [1 update_at]; 
            state.last_R = eye(C);
        end
        Xcov = reshape(Xcov(:,update_at),C,C,[]);
        if usegpu
            Xcov = gather(Xcov); end
        % do the reconstruction in intervals of length stepsize (or shorter if at the end of a chunk)
        last_n = 0;
        for j=1:length(update_at)
            % do a PCA to find potential artifact components
            [V,D] = eig(Xcov(:,:,j));
            [D,order] = sort(reshape(diag(D),1,C)); V = V(:,order);
            % determine which components to keep (variance below directional threshold or not admissible for rejection)
            keep = D<sum((T*V).^2) | (1:C)<(C-maxdims);
            trivial = all(keep);
            % update the reconstruction matrix R (reconstruct artifact components using the mixing matrix)
            if ~trivial
                R = real(M*pinv(bsxfun(@times,keep',V'*M))*V');
            else
                R = eye(C);
            end
            % apply the reconstruction to intermediate samples (using raised-cosine blending)
            n = update_at(j);
            if ~trivial || ~state.last_trivial
                subrange = range((last_n+1):n);
                blend = (1-cos(pi*(1:(n-last_n))/(n-last_n)))/2;
                data(:,subrange) = bsxfun(@times,blend,R*data(:,subrange)) + bsxfun(@times,1-blend,state.last_R*data(:,subrange));
            end
            [last_n,state.last_R,state.last_trivial] = deal(n,R,trivial);
        end
    end
    if splits > 1
        fprintf('.'); end
end
if splits > 1
    fprintf('\n'); end

% carry the look-ahead portion of the data over to the state (for successive calls)
state.carry = [state.carry data(:,(end-P+1):end)];
state.carry = state.carry(:,(end-P+1):end);

% finalize outputs
outdata = data(:,1:(end-P));
if usegpu
    state.iir = gather(state.iir);
    state.cov = gather(state.cov);
end
outstate = state;



function [X,Zf] = moving_average(N,X,Zi)
% Run a moving-average filter along the second dimension of the data.
% [X,Zf] = moving_average(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Channels x #Samples]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2012-01-10

if nargin <= 2 || isempty(Zi)
    Zi = zeros(size(X,1),N); end

% pre-pend initial state & get dimensions
Y = [Zi X]; M = size(Y,2);
% get alternating index vector (for additions & subtractions)
I = [1:M-N; 1+N:M];
% get sign vector (also alternating, and includes the scaling)
S = [-ones(1,M-N); ones(1,M-N)]/N;
% run moving average
X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% read out result
X = X(:,2:2:end);

if nargout > 1
    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end



function result = hlp_memfree
% Get the amount of free physical memory, in bytes
result = java.lang.management.ManagementFactory.getOperatingSystemMXBean().getFreePhysicalMemorySize();
